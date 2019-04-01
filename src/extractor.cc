#include <iostream>
#include <string>
#include <algorithm>
#include <unordered_map>
#include "extractor.h"


using namespace std;

/***************************************************************/
extractor::extractor( string filename, int min_sc, int max_support, int max_fragment_size, double clip_ratio, bool use_indel, bool print_stats): 
	min_sc(min_sc), max_support(max_support), max_fragment_size(max_fragment_size), clip_ratio(clip_ratio), use_indel(use_indel), PRINT_STATS(print_stats)
{

	int min_length = -1;
	FILE *fi = fopen(filename.c_str(), "rb");

	char magic[2];
	fread(magic, 1, 2, fi);
	fclose(fi);

	int ftype = 1; // 0 for BAM and 1 for SAM
	if (magic[0] == char(0x1f) && magic[1] == char(0x8b)) 
		ftype = 0;

	if ( !ftype )
		parser = new BAMParser(filename);
	else
		parser = new SAMParser(filename);

	string comment = parser->readComment();
}

extractor::~extractor(){

	delete parser;
}


/****************************************************************/
int extractor::parse_sc( const char *cigar,  short int &match_l, short int &read_l )
{	
	short int tmp = 0;
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
vector<extractor::breakpoint> extractor::extract_bp(string& cigar, short int& mapped, int sc_loc, bool use_indel){	

	short int len = 0;
	vector<breakpoint> bps;
	mapped = 0;

	for(char c : cigar){

		if(isdigit(c)){
			len *= 10;
			len +=  int(c - '0');
		}
		else{
			if(c == 'M'){
				sc_loc += len;
				mapped += len;
			}
			else if(c != 'D' ){ 
				if(use_indel && c == 'I'){
					bps.push_back({sc_loc,len,BOTH});
				}
				else if(mapped){
					bps.push_back({sc_loc,len,LEFT});
				}
				else{
					bps.push_back({sc_loc,len,RIGHT});
				}
			}
			else{
				if(use_indel)bps.push_back({sc_loc,len,DLEFT});
				sc_loc += len;
				if(use_indel)bps.push_back({sc_loc,len,DRIGHT});
			}

			len = 0;
		}
	}

	return bps;
}

/***************************************************************/
int extractor::dump_oea( const Record &rc, read &tmp, vector<breakpoint> &bps, double clip_ratio )
{
	unordered_map<string, Record>::iterator it;

	it = map_oea.find( rc.getReadName() );
	string cigar, seq = "";
	int flag   = 0, 
		u_flag = 0,
		reversed = 0,
		sc_loc = 0; // if the unmapped end is reversed or not
	short int mapped = 0;
	bool clipped = false;
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

		cigar = string(rc.getCigar());
		sc_loc = rc.getLocation();
		bps = extract_bp(cigar, mapped, sc_loc, false);

		if(bps.empty()){
			for(int i = max_fragment_size/2; i <= max_fragment_size; i++){
				bps.push_back({sc_loc-i, mapped, BOTH});
			}
			for(int i = max_fragment_size/2; i <= max_fragment_size; i++){
				bps.push_back({sc_loc+mapped+i, mapped, BOTH});
			}
		}
			
		//anchor_pos = rc.getLocation(); 
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
int extractor::dump_mapping( const Record &rc, read &tmp, vector<breakpoint> &bps, double clip_ratio )
{
	int flag = 0, u_flag = 0, reversed = 0, sc_loc = 0;
	int flag2 = 0, u_flag2 = 0, reversed2 = 0, sc_loc2 = 0;
	string seq, cigar;
	string seq2, cigar2;

	unordered_map<string, Record >::iterator it = map_read.find( rc.getReadName() );
	if ( it == map_read.end() ) 
	{
		map_read[rc.getReadName()] = rc;
	}
	else
	{	
		Record rc2 = it->second; // with smaller pos
		int part_flag = 0;
		int flag_1 = 0, flag_2 = 0;//, t_flag = 0;
		short r1 = 0, m1 = 0, r2 = 0, m2 = 0; 
		int mate_flag = 0; // 0 for using rc in parition, 1 for using rc2
		int len = 0;
		short int mapped = 0;
		bool clipped = false;

		parse_sc( rc.getCigar(), m1, r1 );
		if ( 0 == r1) { r1 = (short) strlen( rc.getSequence() ); }
		parse_sc( rc2.getCigar(), m2, r2 );
		if ( 0 == r2) { r2 = (short) strlen( rc2.getSequence() ); }

		if ( 0 < r1 && 0 < r2 )
		{
			if ( ( clip_ratio > ( m1 + m2 )*1.0/( r1 + r2 ) ) )
			{  part_flag = 1; }
		}

		if ( part_flag ) // default: use rc2 as anchor and 
		{

			mate_flag = 0;
			reversed = ((rc.getMappingFlag()  & 0x10) == 0x10);
			flag = rc2.getMappingFlag();
			u_flag = rc.getMappingFlag();
			cigar = string(rc.getCigar());
			sc_loc = rc.getLocation();
			//bps = extract_bp(cigar, mapped, sc_loc);
			
			if ( r1 -m1  <  r2 - m2 )
			{
				mate_flag = 1;
				reversed = ((rc2.getMappingFlag()  & 0x10) == 0x10);
				flag = rc.getMappingFlag();
				u_flag = rc2.getMappingFlag();
				cigar = string(rc2.getCigar());
				sc_loc = rc2.getLocation();
				//if(mapped >= 30)bps = extract_bp(cigar, mapped, sc_loc);
			}
			else{
				//cigar = string(rc2.getCigar());
				//sc_loc = rc2.getLocation();
				//extract_bp(cigar, mapped, sc_loc);
				//if(mapped < 30)bps.clear();
			}	

			bps = extract_bp(cigar, mapped, sc_loc, use_indel);

			if (flag & 0x10)
			{ 	//mate has to be positive

				if ( mate_flag ){
					seq = ( reversed ) ? reverse_complement (rc2.getSequence()) : rc2.getSequence();
				}
				else{
					seq = ( reversed ) ? reverse_complement (rc.getSequence()) : rc.getSequence();
				}
				tmp.name = string(rc.getReadName()) + "+";
				tmp.seq = seq;
			}
			else
			{ 
				if ( mate_flag ){
					seq = ( !reversed ) ? reverse_complement (rc2.getSequence()) : rc2.getSequence();
				}
				else{
					seq = ( !reversed ) ? reverse_complement (rc.getSequence()) : rc.getSequence();
				}
				tmp.name = string(rc.getReadName()) + "-";
				tmp.seq = seq;
			}
		}
		else if( 0x2 != ( 0x2 & rc.getMappingFlag() ) ){

			//first mate
			reversed2 = ((rc2.getMappingFlag()  & 0x10) == 0x10);
			flag2 = rc.getMappingFlag();
			u_flag2 = rc2.getMappingFlag();
			cigar2 = string(rc2.getCigar());
			sc_loc2 = rc2.getLocation();
			r2 = (short) strlen( rc2.getSequence() );

			//second mate
			reversed = ((rc.getMappingFlag()  & 0x10) == 0x10);
			flag = rc2.getMappingFlag();
			u_flag = rc.getMappingFlag();
			cigar = string(rc.getCigar());
			sc_loc = rc.getLocation();
			r1 = (short) strlen( rc.getSequence() );

			int diff = (sc_loc-sc_loc2);

			if(reversed == reversed2){
dis_count1++;
				//use first mate
				if(reversed){
					seq = reverse_complement (rc2.getSequence());
					
					for(int i = 0; i < max_fragment_size; i++){
						bps.push_back({sc_loc-i, r1, BOTH});
					}
					
					tmp.name = string(rc2.getReadName()) + "+";
					tmp.seq = seq;
				}
				//use second mate
				else{
					seq = reverse_complement (rc.getSequence());

					for(int i = 0; i < max_fragment_size; i++){
						bps.push_back({sc_loc2+r2+i, r2, BOTH});
					}
					
					tmp.name = string(rc.getReadName()) + "-";
					tmp.seq = seq;
				}

			}
			else if(diff > min_sc && reversed2){
dis_count2++;
				//use both TODO
				seq = rc.getSequence();

				for(int i = 0; i < max_fragment_size; i++){
					bps.push_back({sc_loc2+r2+i, r2, BOTH});
				}
				
				tmp.name = string(rc.getReadName()) + "-";
				tmp.seq = seq;
			}
			else if(diff > min_sc && reversed){
dis_count3++;
			}
		}
		
		map_read.erase( it );
	}
	return (int) map_read.size();
}

/****************************************************************/
// return the entry obtained in the final string
bool extractor::dump_supply( const string& readname, const int flag, read &tmp)
{
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
					tmp.name = readname + "+";
					tmp.seq = it->second.second;
				}
				else if ( !se_flag )
				{
					tmp.name = readname + "-";
					tmp.seq = reverse_complement( it->second.first );
				}
				else
				{
					tmp.name = readname + "-";
					tmp.seq = reverse_complement( it->second.first );
				}
			}
			else // second mate being soft-clipped
			{
				if ( !se_flag ) // forced to place the other mate
				{
					tmp.name = readname + "+";
					tmp.seq = it->second.first;
				}
				else if ( !fr_flag )
				{
					tmp.name = readname + "-";
					tmp.seq = reverse_complement( it->second.second );
				}
				else
				{
					tmp.name = readname + "-";
					tmp.seq = reverse_complement( it->second.second );
				}
			}
		}
		else // clipped reads are forward
		{
			if ( flag & 0x40 ) // first mate being soft-clipped 
			{
				if ( !fr_flag ) // forced to place the other mate
				{
					tmp.name = readname + "-";
					tmp.seq = reverse_complement( it->second.second );
				}
				else if ( !se_flag )
				{
					tmp.name = readname + "+";
					tmp.seq = it->second.first;
				}
				else
				{
					tmp.name = readname + "+";
					tmp.seq = it->second.first;
				}
			}
			else // second mate being clipped
			{
				if ( !se_flag ) // forced to place the other mate
				{
					tmp.name = readname + "-";
					tmp.seq = reverse_complement( it->second.first );
				}
				else if ( !fr_flag )
				{
					tmp.name = readname + "+";
					tmp.seq = it->second.second;
				}
				else
				{
					tmp.name = readname + "+";
					tmp.seq = it->second.second;
				}
			}
		}
		return true;
	}

	return false;
}

/****************************************************************/
void extractor::extract_reads()
{
	vector<breakpoint> bps;
	string readname, cigar;
	read cur_read;
	string read_seq;
	int orphan_flag, oea_flag, chimera_flag;
	int len;
	int sc_loc = 0, start = 0;
	int read_len = 9999999;
	uint32_t flag;
	uint32_t pos;
	uint32_t num_read  = 0;
	char ref[50];
	bool is_supple;
	bool supple_found;
	bool has_supple;
	bool add_read;
	bool cluster_found;
	bool clipped;
	sorted_soft_clips.clear();
	indexed_soft_clips.clear();
	//position_coverage = vector<unsigned short>(300000000, 0);

	while(1)
	{
		if ( parser->hasNext() )
		{
			const Record &rc = parser->next();
			is_supple = false;
			supple_found = false;
			has_supple = false;
			cluster_found = false;
			add_read   = false;
			read_seq = (string)rc.getSequence();
			flag     = rc.getMappingFlag();
			bps.clear();
			sc_loc = 0;

			// for(int i = 0;  i < read_seq.size(); i++){
			// 	position_coverage[rc.getLocation()+i]+=1;
			// }
			if ( MAX_LENGTH < read_seq.size() )
			{
				cerr << "Warning: " << rc.getReadName() <<  " exceeds Maximum Read Length " << MAX_LENGTH << " bp" << endl;
			}

			if ( flag < 256 ) 
			{
				orphan_flag  = (  (flag & 0xc) == 0xc); 
				oea_flag     = ( ((flag & 0xc) == 0x4) || ((flag & 0xc) == 0x8) );
				chimera_flag = ( (0 == (flag & 0xc)  ) && strncmp("=", rc.getPairChromosome(), 1) );
				
				if ( !orphan_flag )
				{
					has_supple = rc.hasSuppleMapping();
					readname = string(rc.getReadName());

					if ( has_supple ){

						auto it    = supply_dict.find(readname);
						if ( it != supply_dict.end() ){
							if ( 0x40 == (flag&0x40) ){
								it->second.first  =  ( 0x10 == (flag&0x10)) ? reverse_complement( read_seq ) : read_seq ;
							}
							else{
								it->second.second  =  ( 0x10 == (flag&0x10)) ? reverse_complement( read_seq ) : read_seq ;
							}
						}
						else{
							if ( 0x40 == (flag&0x40) ){
								supply_dict[readname] = { ( 0x10 == (flag&0x10)) ? reverse_complement( read_seq ) : read_seq, "" }; 
							}
							else{
								supply_dict[readname] = { "", ( 0x10 == (flag&0x10)) ? reverse_complement( read_seq ) : read_seq } ;
							}
						}
					}
				}
				else{
					//orphan_clust.reads.push_back({string(rc.getReadName()), string(rc.getSequence())});
				}

				if ( !orphan_flag && !chimera_flag ){				

					if ( oea_flag ){
oea_count++;
						dump_oea( rc, cur_read, bps, clip_ratio);
					}
					else{
						dump_mapping( rc, cur_read, bps, clip_ratio );			
					}	

					if(!bps.empty()){
						add_read = true;
					}
				}	
			}
			else if ( ( 0x800 == (flag & 0x800) ))
			{		
				// supply_dict does not always include both mate, so we need to check to prevent including empty reads				
				sc_loc   = rc.getLocation();
				cigar    = string(rc.getCigar());
				is_supple = true;
				short int mapped = 0;

				bps = extract_bp(cigar, mapped, sc_loc, use_indel);

				if(!bps.empty()){//} || mapped < 30){
					// insert the hard-clipped mate itself
					supple_found = dump_supply( string(rc.getReadName()), flag, cur_read);
					add_read = true;

					if(!supple_found){
						cur_read.seq =  ( 0x10 == (flag&0x10)) ? reverse_complement( read_seq ) : read_seq;
					}
				}
			}
			
			if ( add_read )
			{
				if(num_read == 0){
					strncpy( ref,  rc.getChromosome(), 50);
					cur_ref = string(ref);
				}
				num_read++;

				if(!start)start = rc.getLocation();

				for(breakpoint& bp : bps){

					if(bp.len >= min_sc){
						sorted_soft_clips.insert({string(rc.getReadName()), cur_read.seq, bp, flag, (is_supple && !supple_found)});
					}
				}

				if(strncmp(ref, rc.getChromosome(), 50 )){
					start = 0;
					num_read = 0;
					if(!sorted_soft_clips.empty()){
						index = sorted_soft_clips.size();
						indexed_soft_clips.reserve(index);
						for(const sortable_read& read : sorted_soft_clips){
							indexed_soft_clips.push_back(read);
						}
						cerr << "chr" << ref << " done" << endl;
						break;
					}
				}
			}

			parser->readNextDiscordant();

		}
		else{
			if(!sorted_soft_clips.empty()){
				index = sorted_soft_clips.size();
				indexed_soft_clips.reserve(index);
				for(const sortable_read& read : sorted_soft_clips){
					indexed_soft_clips.push_back(read);
				}
				cerr << "chr" << ref << " done" << endl;
				cerr << "Processing " << supple_clust.size() << " supplementary clusters..." << endl;
			}
			break;
		}
	}
}

extractor::cluster& extractor::get_next_cluster(int uncertainty, int min_support, bool heuristic) 
{

	if(!supple_clust.empty()){
		supple_clust.pop_back();
	}

	if(index > 0){
		int pos_count, ldel, both, rdel;
		int c_start = 0, c_type;
		extractor::cluster empty_cluster;
		supple_clust.push_back(empty_cluster);

		for(int i = indexed_soft_clips.size()-index; i < indexed_soft_clips.size(); i++){
//if(cur_ref == "9" && i >= 5566 && i <= 6000)cerr << i << "\t" << index << "\t" << indexed_soft_clips[i].bp.sc_loc << endl;
			if(i != indexed_soft_clips.size()-index)index = indexed_soft_clips.size()-i; //probably don't need this.

			sortable_read& sc_read = indexed_soft_clips[i];

			if(cur_pos != sc_read.bp.sc_loc){

				vector<int> counts = vector<int>(5,0);
				counts[sc_read.bp.pos] = 1;
				pos_count = 1;

				while(i+pos_count < indexed_soft_clips.size() && sc_read.bp.sc_loc == indexed_soft_clips[i+pos_count].bp.sc_loc){

					counts[indexed_soft_clips[i+pos_count].bp.pos]++;
					pos_count++;
				}

				ldel = (counts[DLEFT] + counts[LEFT]);
				both = (counts[LEFT] + counts[BOTH] + counts[RIGHT]);
				rdel = (counts[RIGHT] + counts[DRIGHT]);

//if(cur_ref == "7" && sc_read.bp.sc_loc == 55174772)	cerr << counts[DLEFT] << " " << counts[LEFT]  << " " << counts[BOTH] << " " << counts[RIGHT] << " " << counts[DRIGHT] << " for " << i << endl;  95479982
//if(cur_ref == "8" && sc_read.bp.sc_loc >= 42147802 && sc_read.bp.sc_loc <= 42147902)	cerr << counts[DLEFT] << " " << counts[LEFT]  << " " << counts[BOTH] << " " << counts[RIGHT] << " " << counts[DRIGHT] << " for " << i << " @ " << sc_read.bp.sc_loc << endl;


				if(max(max(ldel, both), rdel) < min_support){
					i += (pos_count-1);
					index -= (pos_count+1);
					continue;
				}
				else{
					if(ldel > rdel && ldel > both){
						i--;
						skip_pos = i+ldel;
						skip_count = (rdel+counts[BOTH]);
						cur_type = 0;
					}
					else if(rdel > both){
						i += (ldel+counts[BOTH]-1);
						index -= (ldel+counts[BOTH]);
						cur_type = 2;
					}
					else{
						i += (counts[DLEFT]-1);
						index -= (counts[DLEFT]);
						skip_pos = i+(both+counts[DLEFT]);
						skip_count = counts[DRIGHT];
						cur_type = 1;
					}

					if(skip_count == 0)skip_pos = -1;
					cur_pos = sc_read.bp.sc_loc;
				}
			}
			else if(i == skip_pos){
				i += skip_count;
				index -= (skip_count+1);
				skip_pos = -1;
				skip_count = 0;
				continue;
			}

//if(cur_ref == "7" && sc_read.bp.sc_loc == 55174772) 	cerr << i << " - " << index << "/" << sorted_soft_clips.size() << " ==================== "  << cur_pos << " " << c_start <<  " " << skip_pos << "  " << skip_count << " " << sc_read.seq << endl;
//if(cur_ref == "9" && sc_read.bp.sc_loc == 95479982)cerr << i << " - " << index << "/" << sorted_soft_clips.size() << " ==================== "  << cur_pos << " " << c_start <<  " " << skip_pos << "  " << skip_count << " " << sc_read.seq << endl;
			
			if( !c_start ){
				c_start = sc_read.bp.sc_loc;
				c_type = cur_type;
			}

			if((!heuristic && sc_read.bp.sc_loc != c_start) || (heuristic && uncertainty <= abs(sc_read.bp.sc_loc - c_start))){

				supple_clust.back().start = c_start;
				supple_clust.back().end = c_start;
				supple_clust.back().ref = cur_ref;
				supple_clust.back().total_coverage = 0;//position_coverage[c_start];

				if(!heuristic && !local_reads.empty()){
					sortable_read cur_read = local_reads.front();

					while(uncertainty <= abs(c_start - cur_read.bp.sc_loc)){
						local_reads.pop_front();
						if(local_reads.empty())break;
						cur_read = local_reads.front();
					}

					if(!local_reads.empty()){
						for(auto& local_read : local_reads){

							if((c_start - local_read.bp.sc_loc) == 0)break;

							if(local_read.bp.pos == c_type){
								if(local_read.supple){
									supple_clust.back().sa_reads.push_back((sa_read){local_read.name, local_read.flag});
								}
								else{
									supple_clust.back().reads.push_back({local_read.name, local_read.seq});
								}
							}
						}
					}
				}

				if(supple_clust.back().sa_reads.size() > 0){
					extractor::cluster empty_cluster;
					supple_clust.push_back(empty_cluster);
				}
				else{
					return supple_clust.back();
				}	

				c_start = 0;
			}
			else{
				if(cur_pos != sc_read.bp.sc_loc)i++;
				index--;

				if(!heuristic){
					local_reads.push_back(sc_read);
					local_reads.back().bp.pos = cur_type;
				}

				if(sc_read.supple){
					supple_clust.back().sa_reads.push_back((sa_read){sc_read.name, sc_read.flag});
				}
				else{
					supple_clust.back().reads.push_back({sc_read.name, sc_read.seq});
				}
			}
		}

		supple_clust.back().start = c_start;
		supple_clust.back().end = c_start;
		supple_clust.back().ref = cur_ref;
		supple_clust.back().total_coverage = 0;//position_coverage[c_start];

		if(!heuristic && !local_reads.empty()){
			sortable_read cur_read = local_reads.front();

			while(uncertainty <= abs(c_start - cur_read.bp.sc_loc)){
				local_reads.pop_front();
				if(local_reads.empty())break;
				cur_read = local_reads.front();
			}

			if(!local_reads.empty()){
				for(auto& local_read : local_reads){
					
					if((c_start - local_read.bp.sc_loc) == 0)break;

					if(local_read.bp.pos == c_type){
						if(local_read.supple){
							supple_clust.back().sa_reads.push_back((sa_read){local_read.name, local_read.flag});
						}
						else{
							supple_clust.back().reads.push_back({local_read.name, local_read.seq});
						}
					}
				}
			}
		}

		if(supple_clust.back().sa_reads.size() > 0){
			extractor::cluster empty_cluster;
			supple_clust.push_back(empty_cluster);
			return get_next_cluster(uncertainty, min_support, heuristic); 
		}
		else{
			return supple_clust.back();
		}	

	}
	else if(parser->hasNext()){

		extractor::cluster empty_cluster;
		supple_clust.push_back(empty_cluster);

		extract_reads();

		return get_next_cluster(uncertainty, min_support, heuristic); 
	}	
	else{

		bool add_read;
		read cur_read;

if(PRINT_STATS && supple_clust.size() == 1){
cerr << "Discordant (INV): " << dis_count1 << endl;
cerr << "Discordant (DUP): " << dis_count2 << endl;
cerr << "Discordant (INS): " << dis_count3 << endl;
cerr << "OEA: " << oea_count << endl;
}

		for(auto& read : supple_clust.back().sa_reads){

			add_read = dump_supply( read.name, read.flag, cur_read);

			if (add_read){
				supple_clust.back().reads.push_back({cur_read.name, cur_read.seq});
			}
			else{
				ERROR("Supplemental mappings %s are missing in SAM/BAM, file might be invalid\n", read.name.c_str() );
			}
		}
		
		return supple_clust.back();
	}

}


bool extractor::has_next_cluster(){

	return (parser->hasNext() || (supple_clust.size() > 1));
}

void extractor::clear_maps(){

	supply_dict.clear();
	map_oea.clear();
	map_read.clear();
	map_orphan.clear();

}
