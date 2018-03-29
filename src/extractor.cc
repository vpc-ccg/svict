#include <iostream>
#include <string>
#include <algorithm>
#include <unordered_map>
#include "extractor.h"


using namespace std;

/***************************************************************/
extractor::extractor( string filename, int min_sc, int max_dist, int max_num_read, double clip_ratio ):
	min_sc(min_sc), max_dist(max_dist), max_num_read(max_num_read), clip_ratio(clip_ratio)
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

/****************************************************************/
vector<pair<int, int>> extractor::extract_bp(string& cigar, int& mapped, int sc_loc, bool use_indel){	

	int len = 0;
	vector<pair<int, int>> bps;

	for(char c : cigar){

		if(isdigit(c)){
			len *= 10;
			len +=  int(c - '0');
		}
		else{
			if(c == 'S' || c == 'H' || (use_indel && c == 'I')){
				bps.push_back({sc_loc,len});
			}
			else if(c == 'M'){
				sc_loc += len;
				mapped += len;
			}
			else if(c == 'D' ){ 
				if(use_indel)bps.push_back({sc_loc,len});
				sc_loc += len;
			}
			len = 0;
		}
	}

	return bps;
}

/****************************************************************/
bool extractor::has_supply_mapping( const char *attr )
{
       return ('S' == *attr );
}

/***************************************************************/
int extractor::dump_oea( const Record &rc, read &tmp, vector<pair<int, int>> &bps, double clip_ratio )
{
	unordered_map<string, Record>::iterator it;

	it = map_oea.find( rc.getReadName() );
	string cigar, seq = "";
	int flag   = 0, 
		u_flag = 0,
		reversed = 0,
		sc_loc = 0; // if the unmapped end is reversed or not
	int mapped = 0;
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
			for(int i = INSERT_SIZE/2; i <= INSERT_SIZE; i++){
				bps.push_back({sc_loc-i, mapped});
			}
			for(int i = INSERT_SIZE/2; i <= INSERT_SIZE; i++){
				bps.push_back({sc_loc+mapped+i, mapped});
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
int extractor::dump_mapping( const Record &rc, read &tmp, vector<pair<int, int>> &bps, double clip_ratio )
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
		int r1 = 0, m1 = 0, r2 = 0, m2 = 0; 
		int mate_flag = 0; // 0 for using rc in parition, 1 for using rc2
		int len = 0;
		int mapped = 0;
		bool clipped = false;

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
		else{

			//first mate
			reversed2 = ((rc2.getMappingFlag()  & 0x10) == 0x10);
			flag2 = rc.getMappingFlag();
			u_flag2 = rc2.getMappingFlag();
			cigar2 = string(rc2.getCigar());
			sc_loc2 = rc2.getLocation();
			r2 = (int) strlen( rc2.getSequence() );

			//second mate
			reversed = ((rc.getMappingFlag()  & 0x10) == 0x10);
			flag = rc2.getMappingFlag();
			u_flag = rc.getMappingFlag();
			cigar = string(rc.getCigar());
			sc_loc = rc.getLocation();
			r1 = (int) strlen( rc.getSequence() );

			int diff = (sc_loc-sc_loc2);
			
			if(reversed == reversed2){
				//use first mate
				if(reversed){
					seq = reverse_complement (rc2.getSequence());
					
					for(int i = 0; i < diff; i++){
						bps.push_back({sc_loc2+r2+i, r2});
					}
					
					tmp.name = string(rc2.getReadName()) + "+";
					tmp.seq = seq;
				}
				//use second mate
				else{
					seq = reverse_complement (rc.getSequence());

					for(int i = 0; i < diff; i++){
						bps.push_back({sc_loc-i, r1});
					}
					
					tmp.name = string(rc.getReadName()) + "-";
					tmp.seq = seq;
				}

			}
			else if(diff > min_sc && reversed2){

				//use both TODO
				seq = rc.getSequence();

				for(int i = 0; i < diff; i++){
					bps.push_back({sc_loc-i, r1});
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
	vector<pair<int, int>> bps;
	string readname, cigar;
	read cur_read;
	//read last_read;
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
			flag     = rc.getMappingFlag();
			bps.clear();
			sc_loc = 0;

			if ( flag < 256 ) 
			{
				orphan_flag  = (  (flag & 0xc) == 0xc); 
				oea_flag     = ( ((flag & 0xc) == 0x4) || ((flag & 0xc) == 0x8) );
				chimera_flag = ( (0 == (flag & 0xc)  ) && strncmp("=", rc.getPairChromosome(), 1) );
				
				if ( !orphan_flag )
				{
					has_supple = has_supply_mapping( rc.getOptional() );
					readname = string(rc.getReadName());

					if ( has_supple ){

						auto it    = supply_dict.find(readname);
						if ( it != supply_dict.end() ){
							if ( 0x40 == (flag&0x40) ){
								it->second.first  =  ( 0x10 == (flag&0x10)) ? reverse_complement( rc.getSequence() ) : rc.getSequence() ;
							}
							else{
								it->second.second  =  ( 0x10 == (flag&0x10)) ? reverse_complement( rc.getSequence() ) : rc.getSequence() ;
							}
						}
						else{
							if ( 0x40 == (flag&0x40) ){
								supply_dict[readname] = { ( 0x10 == (flag&0x10)) ? reverse_complement( rc.getSequence() ) : rc.getSequence(), "" }; 
							}
							else{
								supply_dict[readname] = { "", ( 0x10 == (flag&0x10)) ? reverse_complement( rc.getSequence() ) : rc.getSequence() } ;
							}
						}
					}
				}
				else{
					//orphan_clust.reads.push_back({string(rc.getReadName()), string(rc.getSequence())});
				}

				if ( !orphan_flag && !chimera_flag ){				

					if ( oea_flag ){
						dump_oea( rc, cur_read, bps, clip_ratio);
					}
					else{
						dump_mapping( rc, cur_read, bps, clip_ratio );			
					}	

					if(!bps.empty()){// && (cur_read.seq != last_read.seq || cur_read.name == last_read.name)){
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
				int mapped = 0;

				bps = extract_bp(cigar, mapped, sc_loc, use_indel);

				if(!bps.empty()){//} || mapped < 30){
					// insert the hard-clipped mate itself
					supple_found = dump_supply( string(rc.getReadName()), flag, cur_read);
					//if(cur_read.seq != last_read.seq || cur_read.name == last_read.name)
						add_read = true;
				}
			}
			
			if ( add_read )
			{
				//last_read = cur_read;

				if(num_read == 0){
					strncpy( ref,  rc.getChromosome(), 50);
					cur_ref = string(ref);
				}
				num_read++;

				if(!start)start = rc.getLocation();

				if(is_supple && !supple_found)cur_read.seq = "";

				for(pair<int,int>& bp : bps){

					sc_loc = bp.first;

					if(bp.second >= min_sc){
						sorted_soft_clips.push_back({string(rc.getReadName()), cur_read.seq, flag, sc_loc});
						break;
					}
				}

				if(strncmp(ref, rc.getChromosome(), 50 )){
					start = 0;
					num_read = 0;
					if(!sorted_soft_clips.empty()){
						index = sorted_soft_clips.size();
						sort(sorted_soft_clips.begin(), sorted_soft_clips.end());
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
				sort(sorted_soft_clips.begin(), sorted_soft_clips.end());
				break;
			}
		}
	}
}



/****************************************************************/

extractor::cluster& extractor::get_next_cluster_heuristic(int uncertainty, int min_support) 
{

	if(!supple_clust.empty()){
		supple_clust.pop_back();
	}

	if(index > 0){
		int pos_count;
		int cur_pos = -1;
		int num_read = 0;
		int c_start = 0;
		extractor::cluster empty_cluster;
		supple_clust.push_back(empty_cluster);

		for(int i = sorted_soft_clips.size()-index; i < sorted_soft_clips.size(); i++){

			sortable_read& sc_read = sorted_soft_clips[i];
			pos_count = 0;

			if(cur_pos != sc_read.sc_loc){
				while(i+pos_count+1 < sorted_soft_clips.size() && sc_read.sc_loc == sorted_soft_clips[i+pos_count+1].sc_loc){
					pos_count++;
				}

				if(pos_count < min_support){
					i += pos_count;
					index -= (pos_count+1);
					continue;
				}
				else{
					cur_pos = sc_read.sc_loc;
				}
			}

			num_read++;
			
			if( !c_start )c_start = sc_read.sc_loc;

			if(uncertainty < abs(sc_read.sc_loc - c_start)){

				supple_clust.back().start = c_start;
				supple_clust.back().end = c_start;
				supple_clust.back().ref = cur_ref;

				if(supple_clust.back().sa_reads.size() > 0){
					extractor::cluster empty_cluster;
					supple_clust.push_back(empty_cluster);
				}
				else{
					return supple_clust.back();
				}	

				c_start = 0;
				num_read = 0;
			}

			if(sc_read.seq.empty()){
				supple_clust.back().sa_reads.push_back((sa_read){sc_read.name, sc_read.flag});
			}
			else{
				supple_clust.back().reads.push_back({sc_read.name, sc_read.seq});
			}

			index--;
		}

		supple_clust.back().start = c_start;
		supple_clust.back().end = c_start;
		supple_clust.back().ref = cur_ref;

		if(supple_clust.back().sa_reads.size() > 0){
			extractor::cluster empty_cluster;
			supple_clust.push_back(empty_cluster);
			return get_next_cluster(uncertainty, min_support); 
		}
		else{
			return supple_clust.back();
		}	

	}
	else if(parser->hasNext()){

		extractor::cluster empty_cluster;
		supple_clust.push_back(empty_cluster);

		extract_reads();

		return get_next_cluster(uncertainty, min_support); 
	}	
	else{

		bool add_read;
		read cur_read;

		for(auto& read : supple_clust.back().sa_reads){

			add_read = dump_supply( read.readname, read.flag, cur_read);

			if (add_read){
				supple_clust.back().reads.push_back({cur_read.name, cur_read.seq});
			}
			else{
				ERROR("Supplemental mappings %s are missing in SAM/BAM, file might be invalid\n", read.readname.c_str() );
			}
		}
		
		return supple_clust.back();
	}

}

extractor::cluster& extractor::get_next_cluster(int uncertainty, int min_support) 
{

	if(!supple_clust.empty()){
		supple_clust.pop_back();
	}

	if(index > 0){
		int pos_count;
		int cur_pos = -1;
		int num_read = 0;
		int c_start = 0;
		extractor::cluster empty_cluster;
		if(!local_reads.empty()){
			sortable_read cur_read = local_reads.front();
			sortable_read& sc_read = sorted_soft_clips[(sorted_soft_clips.size()-index)];

			while(uncertainty < abs(sc_read.sc_loc - cur_read.sc_loc)){
				local_reads.pop_front();
				if(local_reads.empty())break;
				cur_read = local_reads.front();
			}

			if(!local_reads.empty()){
				for(auto& local_read : local_reads){
					if(local_read.seq.empty()){
						empty_cluster.sa_reads.push_back((sa_read){local_read.name, local_read.flag});
					}
					else{
						empty_cluster.reads.push_back({local_read.name, local_read.seq});
					}
				}
			}
		}
		supple_clust.push_back(empty_cluster);
		

		for(int i = sorted_soft_clips.size()-index; i < sorted_soft_clips.size(); i++){

			sortable_read& sc_read = sorted_soft_clips[i];
			pos_count = 0;

			if(cur_pos != sc_read.sc_loc){

				while(i+pos_count+1 < sorted_soft_clips.size() && sc_read.sc_loc == sorted_soft_clips[i+pos_count+1].sc_loc){
					pos_count++;
				}

				if(pos_count < min_support){
					i += pos_count;
					index -= (pos_count+1);
					continue;
				}
				else{
					cur_pos = sc_read.sc_loc;
				}
			}

			num_read++;
			
			if( !c_start )c_start = sc_read.sc_loc;

			if(sc_read.sc_loc != c_start){

				supple_clust.back().start = c_start;
				supple_clust.back().end = c_start;
				supple_clust.back().ref = cur_ref;

				if(supple_clust.back().sa_reads.size() > 0){

					extractor::cluster empty_cluster;

					if(!local_reads.empty()){
						sortable_read cur_read = local_reads.front();

						while(uncertainty < abs(sc_read.sc_loc - cur_read.sc_loc)){
							local_reads.pop_front();
							if(local_reads.empty())break;
							cur_read = local_reads.front();
						}

						if(!local_reads.empty()){
							for(auto& local_read : local_reads){
								if(local_read.seq.empty()){
									empty_cluster.sa_reads.push_back((sa_read){local_read.name, local_read.flag});
								}
								else{
									empty_cluster.reads.push_back({local_read.name, local_read.seq});
								}
							}
						}
					}

					supple_clust.push_back(empty_cluster);
				}
				else{
					return supple_clust.back();
				}	

				c_start = 0;
				num_read = 0;
			}

			local_reads.push_back(sc_read);

			if(sc_read.seq.empty()){
				supple_clust.back().sa_reads.push_back((sa_read){sc_read.name, sc_read.flag});
			}
			else{
				supple_clust.back().reads.push_back({sc_read.name, sc_read.seq});
			}

			index--;
		}

		supple_clust.back().start = c_start;
		supple_clust.back().end = c_start;
		supple_clust.back().ref = cur_ref;

		if(supple_clust.back().sa_reads.size() > 0){
			extractor::cluster empty_cluster;
			supple_clust.push_back(empty_cluster);
			return get_next_cluster(uncertainty, min_support); 
		}
		else{
			return supple_clust.back();
		}	

	}
	else if(parser->hasNext()){

		extractor::cluster empty_cluster;
		supple_clust.push_back(empty_cluster);

		extract_reads();

		return get_next_cluster(uncertainty, min_support); 
	}	
	else{

		bool add_read;
		read cur_read;

		for(auto& read : supple_clust.back().sa_reads){

			add_read = dump_supply( read.readname, read.flag, cur_read);

			if (add_read){
				supple_clust.back().reads.push_back({cur_read.name, cur_read.seq});
			}
			else{
				ERROR("Supplemental mappings %s are missing in SAM/BAM, file might be invalid\n", read.readname.c_str() );
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