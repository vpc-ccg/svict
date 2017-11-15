#include <string>
#include <map>
#include "common.h"
#include "extractor.h"
#include "record.h"
#include "sam_parser.h"
#include "bam_parser.h"

using namespace std;


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

	fwrite(record.c_str(), 1, record.size(), fp);
}
/****************************************************************/
extractor::extractor(string filename, string output_prefix, int ftype, int oea, int orphan) 
{

	FILE *fi = fopen(filename.c_str(), "rb");
	char magic[2];
	fread(magic, 1, 2, fi);
	fclose(fi);

	string extensions[] = {"","fa","fq","sam"};

	Parser *parser;
	if (magic[0] == char(0x1f) && magic[1] == char(0x8b)) 
		parser = new BAMParser(filename);
	else
		parser = new SAMParser(filename);

	string comment = parser->readComment();


	FILE *foea_mapped, *foea_unmapped, *forphan, *fall_int;
	fall_int  = fopen ((output_prefix + "all_interleaved.fastq").c_str(), "w");
	// open file
	if (oea) 
	{
		string foea_mapped_name = output_prefix + "oea.mapped."+extensions[ftype];
		string foea_unmapped_name = output_prefix + "oea.unmapped."+extensions[ftype];
		foea_mapped  = fopen (foea_mapped_name.c_str(), "w");
		foea_unmapped  = fopen (foea_unmapped_name.c_str(), "w");
	}
		
	if (orphan)
	{
		string forphan_name = output_prefix + "orphan."+extensions[ftype];
		forphan = fopen (forphan_name.c_str(), "w");
	}

	
	uint32_t flag;
	while (parser->hasNext())
	{
		const Record &rc = parser->next();
		flag = rc.getMappingFlag();
		if ((flag & 0xD) == 0xD)
		{
			output_record (forphan, ftype, rc);
			output_record(fall_int, ftype, rc);
		}
		else if ((flag & 0x5) == 0x5) 
		{
			output_record (foea_unmapped, ftype, rc);
			output_record(fall_int, ftype, rc);
		}
		else if ((flag & 0x9) == 0x9)
		{
			output_record (foea_mapped, ftype, rc);
			output_record(fall_int, ftype, rc);
		}
		parser->readNext();
	}
	
	delete parser;
	
	// close file
	if (oea)
	{
		fclose(foea_mapped);
		fclose(foea_unmapped);
	}
	if (orphan)
	{
		fclose(forphan);
	}
	fclose(fall_int);
}
/***************************************************************/
int md_length( char *md)
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
/*************************************************************/
extractor::extractor(string filename, string output) 
{

	FILE *fi = fopen(filename.c_str(), "rb");
	char magic[2];
	fread(magic, 1, 2, fi);
	fclose(fi);

	Parser *parser;
	if (magic[0] == char(0x1f) && magic[1] == char(0x8b)) 
		parser = new BAMParser(filename);
	else
		parser = new SAMParser(filename);

	string comment = parser->readComment();


	FILE *fq = fopen (output.c_str(),"w");
	//FILE *fqsam = fopen ((output + ".sam").c_str(),"w");

	char *opt 			= new char[1000000];
	char *MD 			= new char[1000000];
	string fr_rname;
	string fr_seq;
	string fr_qual;
	uint32_t fr_flag;
	uint32_t flag;
	uint32_t tlen;
	int errNum;
	while (parser->hasNext())
	{
		const Record &rc = parser->next();
		flag = rc.getMappingFlag();
		tlen = strlen(rc.getSequence());
		errNum = int(tlen*0.94);
		if((flag & 0x800) != 0x800)
		{
			if((flag & 0x5) == 0x5 || (flag & 0x9)==0x9 || (flag & 0xD)==0xD)
			{
				output_record(fq, 2, rc);
			}
			else
			{
				strcpy(opt,rc.getOptional());
				strtok(opt,"\t");
				MD = strtok(NULL,"\t");
				if(md_length(MD) >= errNum)
				{
					if((flag & 0x40) == 0x40)
					{
						fr_rname = string(rc.getReadName());
						fr_seq = string(rc.getSequence());
						if((flag & 0x10) == 0x10)
						{
							fr_seq = reverse_complement(fr_seq);
						}
						fr_qual = string(rc.getQuality());
					}
				}
				else
				{
					if((flag & 0x40) == 0x40)
					{
						output_record(fq, 2, rc);
						//output_record(fqsam, 3, rc);
						int secondWritten = 0;
						while(secondWritten != 1)
						{
							parser->readNext();
							if(parser->hasNext())
							{
								const Record &rc2 = parser->next();
								flag = rc2.getMappingFlag();
								if((flag & 0x800) != 0x800)
								{
									output_record(fq, 2, rc2);
									//output_record(fqsam, 3, rc2);
									secondWritten = 1;
								}
							}
						}
					}
					else
					{
						fprintf(fq,"@%s/1\n%s\n+\n%s\n", fr_rname.c_str(), fr_seq.c_str(), fr_qual.c_str());
						//fprintf(fqsam,"@%s/1\t%s\t%s\n", fr_rname.c_str(), fr_seq.c_str(), fr_qual.c_str());
						output_record(fq, 2, rc);
						//output_record(fqsam, 3, rc);
					}
				}
			}
		}
		parser->readNext();
	}
	
	delete parser;
	fclose(fq);
	//fclose(fqsam);
}
/****************************************************************/
int parse_sc( const char *cigar, int &match_l, int &read_l )
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
/****************************************************************/
int process_orphan( const Record &rc, map<string,  Record> &map_orphan, FILE *forphan, FILE *f_int, int ftype)
{
	map<string, Record>::iterator it;
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
int process_oea( const Record &rc, map<string, Record> &map_oea, FILE *f_map, FILE *f_unmap, FILE *f_int, int ftype, int &min_length)
{
	map<string, Record>::iterator it;
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
int examine_mapping( const Record &rc, map<string, Record > &map_read, FILE *f_map, FILE *f_unmap, FILE *f_int, int ftype, double clip_ratio, int &min_length  )
{
	
	map<string, Record >::iterator it = map_read.find( rc.getReadName() );
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
						output_record( f_unmap, ftype, rc2);
						if ( (min_length < 0 ) || ( strlen(rc.getSequence()) < min_length ) ){ min_length = strlen(rc.getSequence());}
					}
					else
					{
						output_record( f_map, ftype, rc2);
						output_record( f_unmap, ftype, rc);
						if ( (min_length < 0 ) || ( strlen(rc2.getSequence()) < min_length ) ){ min_length = strlen(rc2.getSequence());}
					}
					
				}
				else
				{
					if ( mate_flag )
					{
						output_record( f_map, ftype, rc2);
						output_record( f_unmap, ftype, rc);
						if ( (min_length < 0 ) || ( strlen(rc2.getSequence()) < min_length ) ){ min_length = strlen(rc2.getSequence());}
					}
					else
					{
						output_record( f_map, ftype, rc);
						output_record( f_unmap, ftype, rc2);
						if ( (min_length < 0 ) || ( strlen(rc.getSequence()) < min_length ) ){ min_length = strlen(rc.getSequence());}
					}
				}
				
				// interleaved files for genotyping
				if ( ( 0x40 == ( 0x40 & rc.getMappingFlag() ) ) )
				{
					output_record( f_int, 2, rc);
					output_record( f_int, 2, rc2);

				}
				else
				{
					output_record( f_int, 2, rc2);
					output_record( f_int, 2, rc);
				}
			}
		}

		map_read.erase( it );
	}
	return (int) map_read.size();
}
// select any reads whose clipping ratio is less than the clip_ratio in the analysis for mrsfast mapping
/****************************************************************/
extractor::extractor(string filename, string output_prefix, int ftype, int oea, int orphan, double clip_ratio = 0.99 ) 
{
	int min_length = -1;
	FILE *fi = fopen(filename.c_str(), "rb");
	char magic[2];
	fread(magic, 1, 2, fi);
	fclose(fi);

	string extensions[] = {"","fa","fq","sam"};

	Parser *parser;
	if (magic[0] == char(0x1f) && magic[1] == char(0x8b)) 
		parser = new BAMParser(filename);
	else
		parser = new SAMParser(filename);

	string comment = parser->readComment();

	FILE *foea_mapped, *foea_unmapped, *forphan, *fall_int, *f_min_length;
	fall_int  = fopen ((output_prefix + "all_interleaved.fastq").c_str(), "w");
	if (oea) 
	{
		string foea_mapped_name = output_prefix + "oea.mapped."+extensions[ftype];
		string foea_unmapped_name = output_prefix + "oea.unmapped."+extensions[ftype];
		foea_mapped  = fopen (foea_mapped_name.c_str(), "w");
		foea_unmapped  = fopen (foea_unmapped_name.c_str(), "w");
	}
		
	if (orphan)
	{
		string forphan_name = output_prefix + "orphan."+extensions[ftype];
		forphan = fopen (forphan_name.c_str(), "w");
	}

	map<string, Record > map_read;
	map<string, Record > map_orphan;
	map<string,  Record > map_oea;
	int max_size = 0, tmp_size = 0;
	int max_c = 0, tmp_c = 0;
	int max_orphan = 0, max_oea = 0;
	int count = 0;
	
	uint32_t flag;
	uint32_t pos, pair_pos;
	int32_t  tlen;
	int orphan_flag, oea_flag, chimera_flag;
	while (parser->hasNext())
	{
		const Record &rc = parser->next();
		Record cur_rc(rc);

		flag     = rc.getMappingFlag();
		pos      = rc.getLocation();
		pair_pos = rc.getPairLocation();
		tlen     = rc.getTemplateLength();

		if ( flag < 256 ) // To-Do: include supplementary split-mapping as potential mapping locations
		{
			orphan_flag  = (  ( rc.getMappingFlag() & 0xc) == 0xc); 
			oea_flag     = ( (( rc.getMappingFlag() & 0xc) == 0x4) || (( rc.getMappingFlag() & 0xc) == 0x8) );
			chimera_flag = ( (0 == (flag & 0xc)  ) && strncmp("=", rc.getPairChromosome(), 1) );
			
			if ( orphan_flag or chimera_flag) 
			{
				if ( orphan )
				{
					tmp_c = process_orphan( rc, map_orphan, forphan, fall_int, ftype );
					if ( tmp_c > max_orphan) {	max_orphan=tmp_c;	}
				}
			}
			else if ( oea_flag )
			{
				if ( oea )
				{
					tmp_c = process_oea( rc, map_oea, foea_mapped, foea_unmapped, fall_int, ftype, min_length );
					if ( tmp_c > max_oea) {	max_oea = tmp_c; }
				}

			}
			// Hints:  BWA can report concordant mappings with next > pos with negative tlen
			else
			{
				tmp_size = examine_mapping( rc, map_read, foea_mapped, foea_unmapped, fall_int, ftype, clip_ratio, min_length);
				if ( tmp_size > max_size)
				{	
					max_size = tmp_size;
				}
				
			}
			count++; if (0 == count%100000){fprintf( stderr, ".");}
		}
		parser->readNext();
	}
	
	delete parser;
	
	ERROR("");
	//ERROR( "\nMax %d %d %d \n", max_size, max_orphan, max_oea);
	//ERROR( "\nFinal %u %u %u \n", map_read.size(), map_orphan.size(), map_oea.size() );
	//fclose(ftest);
	
	f_min_length  = fopen ((output_prefix + "min_length").c_str(), "w");
	fprintf(f_min_length, "%d\n", min_length);
	fclose( f_min_length);
	// close file
	if (oea)
	{
		fclose(foea_mapped);
		fclose(foea_unmapped);
	}
	if (orphan)
	{
		fclose(forphan);
	}
	fclose(fall_int);
}
/***************************************************************/
extractor::~extractor()
{
		
}
