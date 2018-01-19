#include <iostream>
#include <string>
#include <unordered_map>
#include "extractor.h"


using namespace std;

/***************************************************************/
extractor::extractor( string filename, int max_dist, int max_num_read, double clip_ratio, bool both_mates, bool two_pass ):
	max_dist(max_dist), max_num_read(max_num_read), clip_ratio(clip_ratio), both_mates(both_mates), two_pass(two_pass)
{
	int min_length = -1;
	FILE *fi = fopen(filename.c_str(), "rb");

	char magic[2];
	fread(magic, 1, 2, fi);
	fclose(fi);

	int ftype = 1; // 0 for BAM and 1 for SAM
	if (magic[0] == char(0x1f) && magic[1] == char(0x8b)) 
		ftype = 0;

	if ( two_pass )
	{
		scan_supply_mappings( filename, ftype);
		ERROR("\n%lu supplementary mappings are collected\n", supply_dict.size() );
	}

	if ( !ftype )
		parser = new BAMParser(filename);
	else
		parser = new SAMParser(filename);

	string comment = parser->readComment();
}

extractor::~extractor(){

	delete parser;
}

// 4 is mainly for unmapped format used for building partition
/****************************************************************/
inline void output_record(FILE *fp, int ftype, const Record &rc)
{
	string record;
	uint32_t flag = rc.getMappingFlag();
	
	string mate = ((flag & 0x40)==0x40)?"/1":"/2";
	int reversed = ((flag & 0x10) == 0x10);
	string seq = (reversed) ? reverse_complement (rc.getSequence()) : rc.getSequence();
	string qual = (reversed) ? reverse (rc.getQuality()): rc.getQuality();

	if (ftype == 1)
		record = S(">%s%s\n%s\n", rc.getReadName(), mate.c_str(), seq.c_str());
	else if (ftype==2)
		record = S("@%s%s\n%s\n+\n%s\n", rc.getReadName(), mate.c_str(), seq.c_str(), qual.c_str());
	else if (ftype==3)
		record = S("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\n",
				rc.getReadName(),
				rc.getMappingFlag(),
				rc.getChromosome(),
				rc.getLocation(),
				rc.getMappingQuality(),
				rc.getCigar(),
				rc.getPairChromosome(),
				rc.getPairLocation(),
				rc.getTemplateLength(),
				rc.getSequence(),
				rc.getQuality(),
				rc.getOptional()
				);
	else if ( ftype == 4 )
	{
		record = S("@%s %s\n", rc.getReadName(), seq.c_str() );
	}

	fwrite(record.c_str(), 1, record.size(), fp);
}

/***************************************************************/
int extractor::md_length( char *md)
{
	md+=5;
	int length = 0;
	int tmp = 0;
	while( *md )
	{
		if (isdigit(*md))
		{ 
			tmp = 10 * tmp + (*md - '0');
		}
		else
		{
			length += tmp;
			tmp = 0;
		}
		md++;
	}
	if (0 < tmp){length+=tmp;}
	return length;
}
/****************************************************************/
int extractor::parse_sc( const char *cigar, int &match_l, int &read_l )
{	
	int tmp = 0;
	match_l = 0; read_l = 0;

	while( *cigar )
	{
		if (isdigit(*cigar))
		{ tmp = 10 * tmp + (*cigar - '0');}
		else
		{
			if ( 'M' == *cigar )
			{	match_l += tmp;	}
			if ( 'D' != *cigar )
			{	read_l += tmp;	}
			else
			{ match_l -= tmp; }
			tmp = 0;
		}
		cigar++;
	}
	if (0>match_l){match_l=0;}
	return 0;
}
int extractor::parse_sa( const char *attr )
{
       int flag = 0;
       while( *attr )
       {
               if ('S' == *attr )
               {
                       flag = 1;
               }
               else if ( (1 == flag) && ( 'A' == *attr) )
               {
                       flag = 2;
               }
               else if ( (2 == flag) && ( ':' == *attr) )
               {
                       flag = 3;
                       break;
               }
               else
               {
                       flag = 0;
               }
               attr++;
       }
       return flag;
}

bool extractor::has_supply_mapping( const char *attr )
{
       return ('S' == *attr );
}

/****************************************************************/
int extractor::get_endpoint( const uint32_t pos, const uint32_t pair_pos, const int match_l, const int tlen, int &t_s, int &t_e )
{
	t_s = pos;
	t_e = pos + tlen - 1;
	if ( 0 > tlen )
	{
		t_e = pos + match_l - 1;
		t_s = pair_pos;
	}
}
/****************************************************************/
int extractor::process_orphan( const Record &rc, FILE *forphan, FILE *f_int, int ftype)
{
	unordered_map<string, Record>::iterator it;
	it = map_orphan.find( rc.getReadName() );
	if ( it != map_orphan.end() )	
	{
		if ( ( 0x40 == (0x40 & rc.getMappingFlag() )) )
		{
			output_record( forphan, ftype, rc);
			output_record( forphan, ftype, it->second);
			
			output_record( f_int, 2, rc);
			output_record( f_int, 2, it->second);

		}
		else
		{
			output_record( forphan, ftype, it->second);
			output_record( forphan, ftype, rc);
			
			output_record( f_int, 2, it->second);
			output_record( f_int, 2, rc);
		}
		map_orphan.erase(it);
	}
	else
	{
		map_orphan[ rc.getReadName() ] = rc;
	}
	return (int)map_orphan.size();
}
/****************************************************************/
int extractor::process_oea( const Record &rc, FILE *f_map, FILE *f_unmap, FILE *f_int, int ftype, int &min_length)
{
	unordered_map<string, Record>::iterator it;
	it = map_oea.find( rc.getReadName() );
	if ( it != map_oea.end() )	
	{
		if ( (0x4 == ( 0x4 & rc.getMappingFlag() ) ) )
		{
			output_record( f_unmap, ftype, rc );
			output_record( f_map, ftype, it->second );
			if ( (min_length < 0 ) || ( strlen(it->second.getSequence()) < min_length ) ){ min_length = strlen(it->second.getSequence());}
		}
		else
		{
			output_record( f_map, ftype, rc);
			output_record( f_unmap, ftype, it->second);
			if ( (min_length < 0 ) || ( strlen(rc.getSequence() ) < min_length ) ){ min_length = strlen(rc.getSequence());}
		}
		// interleaved files for genotyping
		if ( ( 0x40 == ( 0x40 & rc.getMappingFlag() ) ) )
		{
			output_record( f_int, 2, rc);
			output_record( f_int, 2, it->second);

		}
		else
		{
			output_record( f_int, 2, it->second);
			output_record( f_int, 2, rc);
		}

		map_oea.erase(it);
	}
	else
	{
		map_oea[ rc.getReadName() ] = rc;
	}
	return (int)map_oea.size();
}

/****************************************************************/
int extractor::examine_mapping( const Record &rc, FILE *f_map, FILE *f_unmap, FILE *f_int, int ftype, double clip_ratio, int &min_length  )
{
	
	unordered_map<string, Record >::iterator it = map_read.find( rc.getReadName() );
	if ( it == map_read.end() ) 
	{
		map_read[rc.getReadName() ] = rc;
	}
	else
	{
		Record rc2 = it->second;
		int flag1 =0, flag2 = 0, t_flag = 0;
		int r1 = 0, m1 = 0, r2 = 0, m2 = 0, 
			tmp_r1  = 0, tmp_m1  = 0,
			tmp_r2 = 0, tmp_m2 = 0;
		int mate_flag = 0; // 0 for rc is first mate; 1 for rc being second mate
		string seq1, qua1, seq2, qua2;

		parse_sc( rc.getCigar(), tmp_m1, tmp_r1 );
		if ( 0 == tmp_r1) { tmp_r1 = (int) strlen( rc.getSequence() ); }
		parse_sc( rc2.getCigar(), tmp_m2, tmp_r2 );
		if ( 0 == tmp_r2) { tmp_r2 = (int) strlen( rc2.getSequence() ); }
		if ( 0x40 == (0x40  & rc.getMappingFlag() ))
		{
			mate_flag = 0;
			if ( r1 <= tmp_r1 && m1 <= tmp_m1 )
			{
				flag1 = rc.getMappingFlag();
				r1    = tmp_r1;
				m1    = tmp_m1;
			} 
			if ( r2 <= tmp_r2 && m2 <= tmp_m2 )
			{
				flag2 = rc2.getMappingFlag();
				r2    = tmp_r2;
				m2    = tmp_m2;
			} 
		}
		else if ( 0x40 == (0x40 & rc2.getMappingFlag() ) )
		{
			mate_flag = 1;
			if ( r2 <= tmp_r1 && m2 <= tmp_m1 )
			{
				flag2 = rc.getMappingFlag();
				r2    = tmp_r1;
				m2    = tmp_m1;
			} 
			if ( r1 <= tmp_r2 && m1 <= tmp_m2 )
			{
				flag1 = rc2.getMappingFlag();
				r1    = tmp_r2;
				m1    = tmp_m2;
			} 
		}

		if ( 0 < r1 && 0 < r2 )
		{
			if ( ( 0x2 != ( 0x2 &flag1) ) || ( clip_ratio > ( m1 + m2 )*1.0/( r1 + r2 ) ) )
			{
				if ( m2 > m1 ) // second mate to mapped, first mate to unmapped
				{
					if ( mate_flag )
					{
						output_record( f_map, ftype, rc);
						output_record( f_unmap, 4, rc2);
						if ( (min_length < 0 ) || ( strlen(rc.getSequence()) < min_length ) ){ min_length = strlen(rc.getSequence());}
					}
					else
					{
						output_record( f_map, ftype, rc2);
						output_record( f_unmap, 4, rc);
						if ( (min_length < 0 ) || ( strlen(rc2.getSequence()) < min_length ) ){ min_length = strlen(rc2.getSequence());}
					}
					
				}
				else
				{
					if ( mate_flag )
					{
						output_record( f_map, ftype, rc2);
						output_record( f_unmap, 4, rc);
						if ( (min_length < 0 ) || ( strlen(rc2.getSequence()) < min_length ) ){ min_length = strlen(rc2.getSequence());}
					}
					else
					{
						output_record( f_map, ftype, rc);
						output_record( f_unmap, 4, rc2);
						if ( (min_length < 0 ) || ( strlen(rc.getSequence()) < min_length ) ){ min_length = strlen(rc.getSequence());}
					}
				}
				
				//// interleaved files for genotyping
				//if ( ( 0x40 == ( 0x40 & rc.getMappingFlag() ) ) )
				//{
				//	output_record( f_int, 2, rc);
				//	output_record( f_int, 2, rc2);

				//}
				//else
				//{
				//	output_record( f_int, 2, rc2);
				//	output_record( f_int, 2, rc);
				//}
			}
		}

		map_read.erase( it );
	}
	return (int) map_read.size();
}

/***************************************************************/
int extractor::dump_oea( const Record &rc, read &tmp, int &anchor_pos, bool both_mates )
{
	unordered_map<string, Record>::iterator it;
	it = map_oea.find( rc.getReadName() );
	string seq = "";
	int flag   = 0, 
		u_flag = 0,
		reversed = 0; // if the unmapped end is reversed or not
	anchor_pos    = 0;
	int mate_flag = 0; // 0 for using rc in parition, 1 for using rc2
	if ( it != map_oea.end() )	
	{
		if ( (0x4 == ( 0x4 & rc.getMappingFlag() ) ) )
		{
			reversed = ((rc.getMappingFlag()  & 0x10) == 0x10);
			//seq = (reversed) ? reverse_complement (rc.getSequence()) : rc.getSequence();
			flag = it->second.getMappingFlag();
			u_flag = rc.getMappingFlag();
		}
		else // decide position
		{
			reversed = ((it->second.getMappingFlag()  & 0x10) == 0x10);
			flag = rc.getMappingFlag();
			u_flag = it->second.getMappingFlag();
			mate_flag = 1;
		}
			
		anchor_pos = rc.getLocation(); 
		if (flag & 0x10)
		{ // anchor reversed, mate positive
			if ( mate_flag )
			{
				seq = ( reversed ) ? reverse_complement (it->second.getSequence()) : it->second.getSequence();
			}
			else
			{
				seq = ( reversed ) ? reverse_complement (rc.getSequence()) : rc.getSequence();
			}
			tmp.name = string(rc.getReadName()) + "+";
			tmp.seq = seq;
		}
		else
		{ 
			if ( mate_flag )
			{
				seq = ( !reversed ) ? reverse_complement (it->second.getSequence()) : it->second.getSequence();
			}
			else
			{
				seq = ( !reversed ) ? reverse_complement (rc.getSequence()) : rc.getSequence();
			}
			tmp.name = string(rc.getReadName()) + "-";
			tmp.seq = seq;
		}
		
		map_oea.erase(it);
	}
	else
	{
		map_oea[ rc.getReadName() ] = rc;
	}
	return (int)map_oea.size();
}

// input: a record and a map for all mappings.
// output: read name along with its location 
/****************************************************************/
int extractor::dump_mapping( const Record &rc, read &tmp, int &anchor_pos, double clip_ratio, bool both_mates )
{
	int flag = 0, u_flag = 0, reversed = 0;;
	anchor_pos = 0;

	unordered_map<string, Record >::iterator it = map_read.find( rc.getReadName() );
	if ( it == map_read.end() ) 
	{
		map_read[rc.getReadName() ] = rc;
	}
	else
	{	
		Record rc2 = it->second; // with smaller pos
		int part_flag = 0;
		int flag_1 = 0, flag_2 = 0;//, t_flag = 0;
		int r1 = 0, m1 = 0, r2 = 0, m2 = 0; 
		int mate_flag = 0; // 0 for using rc in parition, 1 for using rc2
		string seq;


		parse_sc( rc.getCigar(), m1, r1 );
		if ( 0 == r1) { r1 = (int) strlen( rc.getSequence() ); }
		parse_sc( rc2.getCigar(), m2, r2 );
		if ( 0 == r2) { r2 = (int) strlen( rc2.getSequence() ); }
		
		if ( 0 < r1 && 0 < r2 )
		{
			if ( ( 0x2 != ( 0x2 & rc.getMappingFlag() ) ) || ( clip_ratio > ( m1 + m2 )*1.0/( r1 + r2 ) ) )
			{  part_flag = 1; }
		}

		if ( part_flag ) // default: use rc2 as anchor and 
		{
			mate_flag = 0;
			reversed = ((rc.getMappingFlag()  & 0x10) == 0x10);
			//seq = (reversed) ? reverse_complement (rc.getSequence()) : rc.getSequence();
			flag = rc2.getMappingFlag();
			u_flag = rc.getMappingFlag();
			anchor_pos = rc2.getLocation();
			
			if ( r1 -m1  <  r2 - m2 )
			{
				mate_flag = 1;
				reversed = ((rc2.getMappingFlag()  & 0x10) == 0x10);
				flag = rc.getMappingFlag();
				u_flag = rc2.getMappingFlag();
				anchor_pos = rc.getLocation();
			}



			if (flag & 0x10)
			{ 	//mate has to be positive
				//if ( mate_flag )
				//{	seq = (reversed) ? reverse_complement (rc2.getSequence()) : rc2.getSequence();	}
				//else
				//{	seq = (reversed) ? reverse_complement (rc.getSequence()) : rc.getSequence();	}
				//tmp = S("%s+ %s %d\n", rc.getReadName(), seq.c_str(), anchor_pos );
				if ( mate_flag )
				{
					seq = ( reversed ) ? reverse_complement (rc2.getSequence()) : rc2.getSequence();
				}
				else
				{
					seq = ( reversed ) ? reverse_complement (rc.getSequence()) : rc.getSequence();
				}
				tmp.name = string(rc.getReadName()) + "+";
				tmp.seq = seq;
			}
			else
			{ 
				//if ( mate_flag)
				//{	seq = (!reversed) ? reverse_complement (rc2.getSequence()) : rc2.getSequence();	}
				//else
				//{	seq = (!reversed) ? reverse_complement (rc.getSequence()) : rc.getSequence();	}
				//tmp = S("%s- %s %d\n", rc.getReadName(), seq.c_str(), anchor_pos);
				if ( mate_flag )
				{
					seq = ( !reversed ) ? reverse_complement (rc2.getSequence()) : rc2.getSequence();
				}
				else
				{
					seq = ( !reversed ) ? reverse_complement (rc.getSequence()) : rc.getSequence();
				}
				tmp.name = string(rc.getReadName()) + "-";
				tmp.seq = seq;
			}
		}
		
		map_read.erase( it );
	}
	return (int) map_read.size();
}

/****************************************************************/
int extractor::scan_supply_mappings( const string filename, int ftype )
{
	Parser *parser;
	if ( !ftype )
		parser = new BAMParser(filename);
	else
		parser = new SAMParser(filename);

	string comment = parser->readComment();
	int flag = 0;
	bool has_supple = false;
	uint32_t count = 0 ;
	while ( parser->hasNext() )
	{
		const Record &rc = parser->next();
		has_supple = false;

		flag     = rc.getMappingFlag();
		if ( flag < 256 )
		{
			has_supple = has_supply_mapping( rc.getOptional() );// parse SA
			//has_supple = parse_sa( rc.getOptional() );
			if ( has_supple )
			{
				auto it    = supply_dict.find( rc.getReadName() );
				if ( it != supply_dict.end() )
				{
					if ( 0x40 == (flag&0x40) )
					{
						supply_dict[rc.getReadName()].first  =  ( 0x10 == (flag&0x10)) ? reverse_complement( rc.getSequence() ) : rc.getSequence() ;
					}
					else
					{
						supply_dict[rc.getReadName()].second  =  ( 0x10 == (flag&0x10)) ? reverse_complement( rc.getSequence() ) : rc.getSequence() ;
					}
				}
				else
				{
					string token = "";
					if ( 0x40 == (flag&0x40) )
					{
						supply_dict[rc.getReadName()] = { ( 0x10 == (flag&0x10)) ? reverse_complement( rc.getSequence() ) : rc.getSequence(), "" }; 

					}
					else
					{
						supply_dict[rc.getReadName()] = { "", ( 0x10 == (flag&0x10)) ? reverse_complement( rc.getSequence() ) : rc.getSequence() } ;

					}
				}

			}
		}
		count++; if (0 == count%100000){fprintf( stderr, ".");}
		parser->readNextSimple();
	}
	delete parser;
	return 0;
}

/****************************************************************/
// return the entry obtained in the final string
int extractor::dump_supply( const char *readname, const int flag, const size_t pos, bool both_mates, read &tmp)
{
	int num = 0 ;
	auto it    = supply_dict.find( readname );
	if ( it != supply_dict.end() )
	{
		int fr_flag = ( 0 < it->second.first.size()  ) ? 1 : 0;
		int se_flag = ( 0 < it->second.second.size() ) ? 1 : 0;
		// if the read indeed contains first and second mate 

		if ( flag & 0x10 ) // clipped being reversed
		{ 	
			if ( flag & 0x40 ) // first mate being soft-clipped 
			{
				if ( !fr_flag ) // forced to place the other mate
				{
					tmp.name = string(readname) + "+";
					tmp.seq = it->second.second;
					num = 1;
				}
				else if ( !se_flag )
				{
					tmp.name = string(readname) + "-";
					tmp.seq = reverse_complement( it->second.first );
					num = 1;
				}
				else
				{
					tmp.name = string(readname) + "-";
					tmp.seq = reverse_complement( it->second.first );
					num = 2;
				}
			}
			else // second mate being soft-clipped
			{
				if ( !se_flag ) // forced to place the other mate
				{
					tmp.name = string(readname) + "+";
					tmp.seq = it->second.first;
					num = 1;
				}
				else if ( !fr_flag )
				{
					tmp.name = string(readname) + "-";
					tmp.seq = reverse_complement( it->second.second );
					num = 1;
				}
				else
				{
					tmp.name = string(readname) + "-";
					tmp.seq = reverse_complement( it->second.second );
					num = 2;
				}
			}
		}
		else // clipped reads are forward
		{
			if ( flag & 0x40 ) // first mate being soft-clipped 
			{
				if ( !fr_flag ) // forced to place the other mate
				{
					tmp.name = string(readname) + "-";
					tmp.seq = reverse_complement( it->second.second );
					num = 1;
				}
				else if ( !se_flag )
				{
					tmp.name = string(readname) + "+";
					tmp.seq = it->second.first;
					num = 1;
				}
				else
				{
					tmp.name = string(readname) + "+";
					tmp.seq = it->second.first;
					num = 2;
				}
			}
			else // second mate being clipped
			{
				if ( !se_flag ) // forced to place the other mate
				{
					tmp.name = string(readname) + "-";
					tmp.seq = reverse_complement( it->second.first );
					num = 1;
				}
				else if ( !fr_flag )
				{
					tmp.name = string(readname) + "+";
					tmp.seq = it->second.second;
					num = 1;
				}
				else
				{
					tmp.name = string(readname) + "+";
					tmp.seq = it->second.second;
					num = 2;
				}
			}
		}
	}
	else
	{
		ERROR("Supply mappings %s are missing in two-pass scan. SAM/BAM files might be invalid\n", readname );
	}

	return num;
}
/****************************************************************/

extractor::cluster extractor::get_next_cluster() 
{

	extractor::cluster next_cluster;
	vector<pair<string, string>> reads;
	
	int max_c = 0, tmp_c = 0;
	int count = 0;
	
	uint32_t flag;
	uint32_t pos, pair_pos;
	int32_t  tlen;
	int orphan_flag, oea_flag, chimera_flag;

	char ref[1000];
	uint32_t start_loc = 0;
	uint32_t p_loc     = 0;
	uint32_t num_read  = 0;
	uint32_t num_mappings  = 0;
	uint32_t base      = 0;
	int      index     = 0;

	read tmp;	
	int s1 = 0, e1 = 0;
	int match_l = 0, read_l = 0;
	int t_s = 0, t_e = 0;
	int t_loc;

	int cluster_id = 1;
	fpos_t cur_pos;
	int p_start = 0, p_end = 0;
	int cluster_flag = 1;
	
	vector< read > vec_read;
	vec_read.reserve(max_num_read);

	int num_supply = 0;
	int fr_flag = 0, se_flag = 0; // existence for first and second mate in supply_dict
	int add_read = 0 ;
	//int reversed = 0, // 1 for rc, 0 for forward
	//	mate = 0; 	  // 1 for first mate, 0 for second	
	//int z = 0;
	while(1){
		if ( parser->hasNext() )
		{

			const Record &rc = parser->next();

			add_read   = 0;

			flag     = rc.getMappingFlag();
			pos      = rc.getLocation();

			if(vec_read.empty())strncpy( ref,  rc.getChromosome(), 1000);

			if ( flag < 256 ) // To-Do: include supplementary split-mapping as potential mapping locations
			{
				orphan_flag  = (  (flag & 0xc) == 0xc); 
				oea_flag     = ( ((flag & 0xc) == 0x4) || ((flag & 0xc) == 0x8) );
				chimera_flag = ( (0 == (flag & 0xc)  ) && strncmp("=", rc.getPairChromosome(), 1) );

				if ( !orphan_flag and !chimera_flag )
				{
					t_loc = 0; 
					if ( oea_flag )
					{
						tmp_c = dump_oea( rc, tmp, t_loc, both_mates);
					}
					else
					{
						tmp_c = dump_mapping( rc, tmp, t_loc, 0.99, both_mates );			
					}	
					
					if ( t_loc )
					{
						add_read = (both_mates) ? 2: 1;

					}

				}
			}
			else if ( two_pass && ( 0x800 == (flag & 0x800) ) )
			{		
				// supply_dict does not always include both mate, so we need to check to prevent including empty reads
				//
				// insert the hard-clipped mate itself
				add_read  = dump_supply( rc.getReadName(), flag, pos, both_mates, tmp);
				if ( add_read )
				{ 
					t_loc = pos;
				}
			}

			if ( add_read )
			{

				vec_read.push_back( tmp );  
				num_read++; 
				num_mappings += add_read;
				p_end = t_loc;
				if ( !p_start ){ p_start = t_loc;}

				if ( strncmp(ref, rc.getChromosome(), 1000 ) || ( max_dist < t_loc - p_start)  || ( max_num_read < num_read) )
				{
					if ( num_read )//&& cluster_flag )
					{

						for (auto &i: vec_read)
							reads.push_back({i.name,i.seq});	

						next_cluster.reads = reads;
						next_cluster.start = p_start;
						next_cluster.end = p_end;
						next_cluster.ref = string(ref);
					}

					p_start     = 0;
					p_end       = 0;
					num_read    = 0;
					num_mappings = 0;
					strncpy( ref,  rc.getChromosome(), 1000);
					parser->readNextSimple();
					break;
				}
			}
			
			parser->readNextSimple();
		}
		else{
			next_cluster.reads = reads; 
			break;
		}
	}

	return next_cluster;

}

bool extractor::has_next_cluster(){

	return parser->hasNext();
}