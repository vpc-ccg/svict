#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include "annotation.h"

int _sample_num = 0;

int ana_flag = 0;

using namespace std;
/**********************************************/
bool comp_isoform_start( const isoform &is_1, const isoform &is_2)
{
	int l1 = (int) is_1.exon.size();
	int l2 = (int) is_2.exon.size();
	//E("%s %s %d %d\n", is_1.id.c_str(), is_2.id.c_str(), l1, l2);
	return ( (is_1.exon[0].e_s < is_2.exon[0].e_s) || ( is_1.exon[l1-1].e_e < is_2.exon[l2-1].e_e));
}

/**********************************************/
bool comp_exon( const Exon &exon1, const Exon &exon2)
{
	return ( exon1.e_s < exon2.e_s );
}

/**********************************************/
bool comp_region( const Region &r1, const Region &r2)
{
	return ( (r1.start < r2.start) || ( (r1.start == r2.start) && (r1.end < r2.end) ) );
}
/**********************************************/
bool comp_gene_in_contig( const gene_data &g1, const gene_data &g2)
{
	return ( (g1.start < g2.start) || ( (g1.start == g2.start) && (g1.end < g2.end) ) );
}

/**********************************************/
uint32_t overlap_l( const uint32_t &s1, const uint32_t &e1, const uint32_t &s2, const uint32_t &e2)
{
	uint32_t l = 0;
	if ( s1 <= s2 && e2 <= e2 )
	{ l = e2 - s2 + 1;	}
	else if ( s2 <= s1 && e1 <= e2 )
	{	l = e1 - s1 + 1;	}
	else if ( s1 <= s2 && s2 <= e1)
	{	l = e1 - s2 + 1;	}
	else if ( s2 <= s1 && s1 <= e2)
	{	l =	e2 - s1 + 1;	}
	return l;
}
/**********************************************/
void adjust_gene( map<string, gene_data> &map_gene, const string g_id, char *gene_name, const isoform &new_iso )
{
	map<string, gene_data>::iterator it;
	it = map_gene.find(g_id);

	uint32_t left = 0, right = 0;
	int limit = (int)new_iso.exon.size();
	left = new_iso.exon[0].e_s;
	right = new_iso.exon[limit-1].e_e;

	if ( it != map_gene.end())
	{
		if ( it->second.start > left )
		{
			it->second.start = left;
		}
		if ( it->second.end < right )
		{
			it->second.end = right;
		}
		if ( 1 == new_iso.cds_iso)
		{
			it->second.cds_gene = 1;
		}
	}
	else
	{
		gene_data tmp_gene;
		tmp_gene.gene_id = g_id;
		tmp_gene.gene_name = string(gene_name);
		tmp_gene.chr = new_iso.ref;
		tmp_gene.strand = new_iso.strand;
		tmp_gene.start = left;
		tmp_gene.end = right;
		tmp_gene.cds_gene = new_iso.cds_iso;
		map_gene[ g_id ] = tmp_gene;
	}	
}
/**********************************************/
// We are interested in either CDS or exons. UTR and start/stop codon can be skipped for now.
int feature_of_interest( char *fea ) 
{
	int flag = ( !strncmp( "C", fea, 1) || !strncmp("e", fea, 1) );
	return flag;
}
/**********************************************/
// Reading Ensembl Version GTF File with the following properties:
// Entries of same isoform are grouped together (including both exons, CDS, and UTR); 
// Entries of the same gene are grouped together;
// ToDo: We do not fully make use of these properties in GRCh37.75 
void ensembl_Reader(const char *gtf_file, map<string, vector <isoform> > &iso_gene_map, map<string, vector<gene_data> > &gene_sorted_map)
{
	bool new_gtf 	= false;
	char *readline 	= (char*)malloc(MAX_LINE_ANNO);
	char *seq 		= (char*)malloc(TOKEN_LENGTH_ANNO);
	char *src	 	= (char*)malloc(TOKEN_LENGTH_ANNO);
	char *fea 		= (char*)malloc(TOKEN_LENGTH_ANNO);
	char *misc 		= (char*)malloc(TOKEN_LENGTH_ANNO);
	char *gid_str 	= (char*)malloc(TOKEN_LENGTH_ANNO);
	char *gn_str 	= (char*)malloc(TOKEN_LENGTH_ANNO);
	char *tid_str 	= (char*)malloc(TOKEN_LENGTH_ANNO);
	char *tn_str 	= (char*)malloc(TOKEN_LENGTH_ANNO);
	char *bio_str 	= (char*)malloc(TOKEN_LENGTH_ANNO);
	
	uint32_t start = 0, end = 0; 
	// start and ending position of exon unit and cds region
	uint32_t exon_s = 0, exon_e = 0, cds_s = 0, cds_e = 0; 
	// uint32_t s1 = 0, e1 = 0, s2 = 0, e2 = 0; // start and ending position of exon unit and cds region
	
	char strand;
	char del;
	char fr_cha;
	int  frame=0, offset=0, of2=0;
	int  count = 0;	
	int  len = 0;
	string gstr;
	char *gene_id 	 = (char*)malloc(TOKEN_LENGTH_ANNO);
	char *gene_name  = (char*)malloc(TOKEN_LENGTH_ANNO);
	char *trans_id 	 = (char*)malloc(TOKEN_LENGTH_ANNO);
	char *trans_name = (char*)malloc(TOKEN_LENGTH_ANNO);
	trans_id[0]='\0';
	CDS cds_token;	
	Exon exon_token;
	isoform iso_token;
	//GTFAttr gtf_attr;

	int cds_gene = 0; // 1 if any transcripts of the gene code protein
	map< string, gene_data > map_gene; // gene id to gene data
	gene_data g_token;

	int num_attr = 0;
	FILE *fp = fopen( gtf_file, "r" );
	while( NULL != fgets( readline, MAX_LINE_ANNO, fp) )
	{
		if ( 0 == strncmp("#", readline, 1) ) {	

			char * token = strtok(readline," ");
			if(0 == strcmp("#!genome-version", token)){
				token = strtok(NULL, " ");
				char ver[3];
				memcpy(ver, &token[4], 2);
				ver[3] = '\0';
				int ver_num = atoi(ver);
				if(ver_num >= 38)new_gtf = true;
			}
			continue; 
		}
		
		sscanf( readline, "%s %s %s %u %u %s %c %c%n",
			seq, src, fea,  &start, &end, misc, &strand, &fr_cha, &offset );


//1	protein_coding	exon	907455	907530	.	+	.	gene_id "ENSG00000187583"; 					  transcript_id "ENST00000379407"; 						   exon_number "9"; gene_name "PLEKHN1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PLEKHN1-004"; transcript_source "ensembl_havana"; tag "CCDS"; ccds_id "CCDS53256"; exon_id "ENSE00001375540";
//1	ensembl			exon	939275	939412	.	+	.	gene_id "ENSG00000187634"; gene_version "11"; transcript_id "ENST00000617307"; transcript_version "4"; exon_number "7"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTHUMG00000040719"; havana_gene_version "10"; transcript_name "SAMD11-203"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "ENSE00001708361"; exon_version "1"; tag "basic"; transcript_support_level "5";

		
		if ( !feature_of_interest( fea ) ) { continue; }

		if(new_gtf){
			// exon and CDS only
			sscanf( readline + offset, 
			"%s %c%[^\"]%s\
			 %s %s\
			 %s %c%[^\"]%s\ 
			 %s %s\
			 %s %s\
			 %s %c%[^\"]%s %n",
			misc, &del, gid_str, misc, 
			misc, misc,    
			misc, &del, tid_str, misc,
			misc, misc,  
			misc, misc,
			misc, &del, gn_str, misc, &of2);
		}
		else{
			// exon and CDS only
			sscanf( readline + offset, 
			"%s %c%[^\"]%s\
			 %s %c%[^\"]%s\ 
			 %s %s\
			 %s %c%[^\"]%s %n",
			misc, &del, gid_str, misc,   
			misc, &del, tid_str, misc,
			misc, misc,
			misc, &del, gn_str, misc, &of2);
		}


		
		// extra string copy slows down the following function
		// num_attr = get_gtf_info( readline + offset, gtf_attr );
		//L(">>GID%s\n", gid_str);L("$$TID%s\n", tid_str);L("^^GN%s\n", gn_str); //L("<<%s\n", tn_str); //L("--%s\n", bio_str); 
		
		// Start A New Transcript
		if ( 0 != strncmp( trans_id, tid_str, TOKEN_LENGTH_ANNO ) )
		{
			if ( '\0' != trans_id[0] )
			{	
				if ( 0 < exon_s + exon_e )	// last region is not a cds??? What r u talking about
				{	
					exon_token.e_s = exon_s;
					exon_token.e_e = exon_e;
					exon_token.c_s = cds_s;
					exon_token.c_e = cds_e;
					//// determine if a token is cds, utr, or intron
					exon_token.cds = 0;
					if ( 0 < cds_s + cds_e )
					{
						exon_token.cds = 2;
						if ( ( exon_s != cds_s) || ( exon_e != cds_e) )
						{	exon_token.cds = 1;	}
					}

					len += (exon_e - exon_s + 1);
					iso_token.exon.push_back(exon_token);	
				}
				//L("Finish Isoform %s with %d exons\n", trans_id, (int)iso_token.exon.size() ); 
				gstr = string(gene_id);				
				sort( iso_token.exon.begin(), iso_token.exon.end(), comp_exon);
				if ( 1 == cds_gene ){ iso_token.cds_iso = 1; }
				iso_gene_map[gstr].push_back( iso_token );

				if ( 0 == iso_token.exon.size()) { E("Warning: No Exons in %s\n", iso_token.id.c_str() ) ;}
				else{ adjust_gene( map_gene, gstr, gene_name, iso_token ); }
			}
			//L("Starting New Isoform %s\t%s\n", tid_str, gid_str);
			strncpy( gene_id,    gid_str, TOKEN_LENGTH_ANNO );
			strncpy( gene_name,  gn_str,  TOKEN_LENGTH_ANNO );
			strncpy( trans_id,   tid_str, TOKEN_LENGTH_ANNO );
			strncpy( trans_name, tn_str,  TOKEN_LENGTH_ANNO );
			iso_token.id      = string(trans_id);
			iso_token.tname   = string(trans_id);
			iso_token.gid     = string(gene_id);
			iso_token.ref     = string(seq);
			iso_token.src     = string(src);
			iso_token.fea     = string(fea);
			iso_token.strand  = strand;
			iso_token.cds_iso = 0;
			iso_token.len     = len;
			iso_token.exon.clear();
			exon_s = 0; exon_e = 0; cds_gene = 0; len = 0;
		}

		// CDS and Exon
		// Note in Ensembl GTF, CDS comes after the corresponding Exon record
		//
		if ( 0 == strncmp(fea, "exon", 4) )
		{	// Previous records
			exon_token.e_s = exon_s;
			exon_token.e_e = exon_e;
			exon_token.c_s = cds_s;
			exon_token.c_e = cds_e;
			//// determine if a token is cds, utr, or intron
			exon_token.cds = 0;
			if ( 0 < cds_s + cds_e )
			{
				exon_token.cds = 2;
				if ( ( exon_s != cds_s) || ( exon_e != cds_e) )
				{	exon_token.cds = 1;	}
			}
			if ( 0 < exon_s + exon_e )	
			{	
				len += (exon_e - exon_s + 1);
				iso_token.exon.push_back(exon_token);
			}

			
			// Reading Current Records
			exon_s = start;
			exon_e = end;
			cds_s = 0;
			cds_e = 0;
		}
		else if ( 0 == strncmp(fea, "CDS", 3) )
		{
			cds_s = start;
			cds_e = end;
			cds_gene = 1; // reading a protein-coding transcript
		}
		
		count++;
		if ( 0 == count%1000000){E(".");}

	}
	//Last Record
	exon_token.e_s = exon_s;
	exon_token.e_e = exon_e;
	exon_token.c_s = cds_s;
	exon_token.c_e = cds_e;
	exon_token.cds = 0;
	if ( 0 < cds_s + cds_e )
	{
		exon_token.cds = 2;
		if ( ( exon_s != cds_s) || ( exon_e != cds_e))
		{	exon_token.cds = 1;	}
	}
	if ( 0 < exon_s + exon_e )	
	{	iso_token.exon.push_back(exon_token);	}
	//L("Finish Isoform %s with %d exon\n", trans_id, (int)iso_token.exon.size() );
	//strncpy( gene_name, gn_str, TOKEN_LENGTH_ANNO);
	
	gstr = string(gene_id);
	sort( iso_token.exon.begin(), iso_token.exon.end(), comp_exon);
	if ( 1 == cds_gene ){ iso_token.cds_iso = 1;}
	iso_gene_map[gstr].push_back(iso_token);
	if ( 0 == iso_token.exon.size()) { E("Warning: No Exons in %s\n", iso_token.id.c_str() ) ;}
	else{adjust_gene( map_gene, gstr, gene_name, iso_token );}
	
	E("\n");
	free(readline);
	free(seq);
	free(src);
	free(fea);
	free(misc);
	free(gid_str);
	free(gn_str); 
	free(tid_str);
	free(tn_str);
	free(bio_str);
	free(gene_id);
	free(gene_name);
	free(trans_id);
	free(trans_name);

	map<string, gene_data >::iterator it;
	for( it = map_gene.begin(); it != map_gene.end(); it++)
	{
		//L("Add_Gene\t%s\t%s\t%s\t%u\t%u\t%d\n", it->second.chr.c_str(), it->second.gene_id.c_str(), it->second.gene_name.c_str(), it->second.start, it->second.end, it->second.cds_gene );
		gene_sorted_map[it->second.chr].push_back(it->second);
	}
	E("Scanning total %d genes\n", (int)map_gene.size() );
	
	map<string, vector<gene_data> >::iterator git;
	for( git = gene_sorted_map.begin(); git != gene_sorted_map.end(); git++ )
	{
		sort( git->second.begin(), git->second.end(), comp_gene_in_contig);
	}
}


/**********************************************/
// Output Gene/Isoform Boundary to outfile for masking purpose
void output_Isoform_Boundary( map<string, vector <isoform> > &iso_gene_map, char *outfile)
{
	//E("ok");
	FILE *out_fp = fopen( outfile, "w" );
	map<string, vector<isoform> >::iterator it;
	for ( it = iso_gene_map.begin(); it != iso_gene_map.end(); it++)
	{
		L("chr %s\n", it->first.c_str());
		int limit = it->second.size();
		L("chr_count %d\n", limit);
		//sort( it->second.begin(), it->second.end(), comp_isoform_start);
		for ( int i = 0; i < limit; i++)
		{
			L("start %d\n", i);
			int iso_l = it->second[i].exon.size();
			L("size %d\n", iso_l);
			
			fprintf(out_fp,"%s\t%u\t%u\tgene_id\"%s\";transcript_id\"%s\";src\"%s\";\n",it->second[i].ref.c_str(), it->second[i].exon[0].e_s,  it->second[i].exon[iso_l -1].e_e, 
				it->second[i].gid.c_str(), it->second[i].id.c_str(), it->second[i].src.c_str());
		}
		
	}

	fclose(out_fp);
	
}
/**********************************************/
void bed_reader( char *in_file, map< string, vector < Region > > &region_map )
{
	char *readline = (char*)malloc(MAX_LINE_ANNO);
	char *ref      = (char*)malloc(TOKEN_LENGTH_ANNO);
	char *src = (char*)malloc(TOKEN_LENGTH_ANNO);
	char *fea = (char*)malloc(TOKEN_LENGTH_ANNO);
	char *misc = (char*)malloc(TOKEN_LENGTH_ANNO);
	char *gid_str = (char*)malloc(TOKEN_LENGTH_ANNO);
	char *gn_str = (char*)malloc(TOKEN_LENGTH_ANNO);
	char *tid_str = (char*)malloc(TOKEN_LENGTH_ANNO);
	char *tn_str = (char*)malloc(TOKEN_LENGTH_ANNO);
	
	int  offset=0;
	uint32_t start =0, end = 0; 
	
	Region t_region;
	map<string, vector <Region> >::iterator it;	
	FILE *fp = fopen( in_file, "r" );
	while( NULL != fgets( readline, MAX_LINE_ANNO, fp) )
	{
		sscanf(readline, "%s %u %u %n", ref, &start, &end, &offset);
		t_region.start = start;
		t_region.end   = end;
		it = region_map.find( ref );
		if ( it == region_map.end() )
		{
			vector<Region> t_vec;
			t_vec.reserve(1024);
			region_map[ ref] = t_vec;
		}
		region_map[ref].push_back( t_region );
	}
	fclose(fp);
	
}
/**********************************************/
void  locate_in_isoform( uint32_t s, uint32_t e, const isoform &iso, vector<uint32_t> &t_vec, bool genomic)
{
	int limit = (int)iso.exon.size();
	uint32_t exon_match = 0, cds_match = 0;
	uint32_t e_s = 0, e_e = 0, c_s = 0, c_e = 0; // exon s/e and cds s/e for the best overlapping record
	uint32_t t_e, t_c;
	uint32_t r_s = 0, r_e = 0;
	uint32_t tr_s = 0, tr_e = 0, exon_len = 0;

	if(!genomic && iso.strand == '-'){
	//	s = iso.len - e + 1;
	//	e = iso.len - s + 1;
	}

//cerr << "start" << endl;
	for (int i = 0 ; i < limit; i ++)
	{
		exon_len = (iso.exon[i].e_e - iso.exon[i].e_s + 1);
		tr_s = tr_e;
		tr_e += exon_len;
		
		if(genomic){
			if ( iso.exon[i].e_e < s )
			{
				r_s += exon_len;
				r_e = r_s;
				continue;
			}
			else if ( e < iso.exon[i].e_s)
			{	continue; }//break;	}
			else
			{
				r_s += (s - iso.exon[i].e_s + 1);
				r_e += e < iso.exon[i].e_e ? (e - iso.exon[i].e_s + 1) : exon_len;
				t_e = (int)overlap_l( s, e, iso.exon[i].e_s, iso.exon[i].e_e);
				t_c = (int)overlap_l( s, e, iso.exon[i].c_s, iso.exon[i].c_e);
				if ( exon_match < t_e )
				{ exon_match = t_e; e_s = iso.exon[i].e_s; e_e = iso.exon[i].e_e;}
				if ( cds_match < t_e )
				{ cds_match = t_c;  c_s = iso.exon[i].c_s; c_e = iso.exon[i].c_e;}
			}
		}
		else{

			if ( tr_e < s )
			{
				continue;
			}
			else if ( e < tr_s)
			{	continue; }//break;	}
			else
			{
				if(iso.strand == '-'){
					r_s = iso.exon[i].e_e - (s - tr_s - 1) ;
					r_e = e < tr_e ? (iso.exon[i].e_e  - (e - tr_s - 1)) : iso.exon[i].e_e;
				}
				else{
					r_s = iso.strand == '-' ? iso.exon[i].e_e - (s - tr_s - 1) : iso.exon[i].e_s + (s - tr_s - 1);
					r_e = e < tr_e ? (iso.exon[i].e_s  + (e - tr_s - 1)) : iso.exon[i].e_e;
				}

				t_e = (int)overlap_l( r_s, r_e, iso.exon[i].e_s, iso.exon[i].e_e);
				t_c = (int)overlap_l( r_s, r_e, iso.exon[i].c_s, iso.exon[i].c_e);
				if ( exon_match < t_e )
				{ exon_match = t_e; e_s = iso.exon[i].e_s; e_e = iso.exon[i].e_e;}
				if ( cds_match < t_e )
				{ cds_match = t_c;  c_s = iso.exon[i].c_s; c_e = iso.exon[i].c_e;}
			}
		}
	}
//cerr << "end" << endl;
	if(genomic && iso.strand == '-'){
		r_s = iso.len - r_e + 1;
		r_e = iso.len - r_s + 1;
	}

	t_vec.clear();
	t_vec.reserve(8);
	t_vec.push_back(exon_match);	t_vec.push_back(cds_match);
	t_vec.push_back(e_s);	t_vec.push_back(e_e);
	t_vec.push_back(c_s);	t_vec.push_back(c_e);
	t_vec.push_back(r_s);	t_vec.push_back(r_e);
}
/**********************************************/
int locate_interval( const string &ref, uint32_t s, uint32_t e, const vector<gene_data> &gene_vector, int pos, const map<string, vector< isoform> > &iso_gene_map,
	string &best_gene, string &best_name, string &best_trans, vector<uint32_t> &vec_best)
{
	map<string, vector<uint32_t>> vec_all;

	return locate_interval(ref, s, e, gene_vector, pos, iso_gene_map, best_gene, best_name, best_trans, vec_best, vec_all);
}
/**********************************************/
int locate_interval( const string &ref, uint32_t s, uint32_t e, const vector<gene_data> &gene_vector, int pos, const map<string, vector< isoform> > &iso_gene_map,
	string &best_gene, string &best_name, string &best_trans, vector<uint32_t> &vec_best, map<string, vector<uint32_t>> &vec_all)
{
	uint32_t max_e = 0, max_c = 0; // maximum overlap between [s,e] with any exons and any CDSs
	uint32_t e_s = 0, e_e = 0, c_s = 0, c_e = 0;
	vector<uint32_t> vec_match;
	best_gene  = "intergenic";
	best_name  = "intergenic";
	best_trans = "intron";


	int flag = 0; // 0 for intergenic, 1 for intronic, 2 for exonic-NC, 3 for UTR, 4 for CDS

	int limit = gene_vector.size();
	map<string, vector<isoform> >::const_iterator iso_it;
	for ( int i = pos; i < limit; i++ )
	{
		if ( gene_vector[i].end < s  )
		{
			pos ++;
			continue;
		}
		else if ( e < gene_vector[i].start )
		{	break;	}
		else
		{
			string gene_id = gene_vector[i].gene_id;
			string gene_name = gene_vector[i].gene_name;
			if (  0 ==  flag){ 
				flag = 1; 
				best_gene = gene_id;
				best_name = gene_name;
			}
			
			iso_it = iso_gene_map.find( gene_id );
			int num_iso = iso_it->second.size();
			for( int j = 0; j < num_iso; j++)
			{
				locate_in_isoform( s, e, iso_it->second[j], vec_match, true);

				vec_all[iso_it->second[j].id].reserve(3);
				vec_all[iso_it->second[j].id].push_back( (uint32_t)flag); //TODO this flag would be incorrect
				vec_all[iso_it->second[j].id].push_back(vec_match[6]);	vec_all[iso_it->second[j].id].push_back(vec_match[7]);

				if (max_e < vec_match[0])
				{
					max_e = vec_match[0];
					e_s   = vec_match[2];	e_e   = vec_match[3];
					best_trans = iso_it->second[j].id;
					if( 2 > flag) {flag = 2 + iso_it->second[j].cds_iso;}
				}
				if (max_c < vec_match[1])
				{
					max_c = vec_match[1];
					c_s   = vec_match[4];	c_e   = vec_match[5];
					best_trans = iso_it->second[j].id;
					flag = 4;
				}
			}
		}
	}

	if(best_gene == "intergenic")best_trans = best_gene;

	vec_best.clear();
	vec_best.reserve(8);
	vec_best.push_back( (uint32_t)flag);
	vec_best.push_back(max_e);	vec_best.push_back(max_c);
	vec_best.push_back(e_s);	vec_best.push_back(e_e);
	vec_best.push_back(c_s);	vec_best.push_back(c_e);

	return pos;
}

/**********************************************/
int locate_interval( const string &ref, uint32_t s, uint32_t e, const string &gene_id, const string &trans_id, const vector<gene_data> &gene_vector, int pos, const map<string, vector< isoform> > &iso_gene_map,
	string &best_gene, string &best_name, string &best_trans, vector<uint32_t> &vec_best, map<string, vector<uint32_t>> &vec_all)
{
	uint32_t max_e = 0, max_c = 0; // maximum overlap between [s,e] with any exons and any CDSs
	uint32_t e_s = 0, e_e = 0, c_s = 0, c_e = 0;
	vector<uint32_t> vec_match;
	best_gene  = "intergenic";
	best_name  = "intergenic";
	best_trans = "intron";


	int flag = 0; // 0 for intergenic, 1 for intronic, 2 for exonic-NC, 3 for UTR, 4 for CDS


	map<string, vector<isoform> >::const_iterator iso_it;

	best_gene = gene_id;

	iso_it = iso_gene_map.find( gene_id );
	int num_iso = iso_it->second.size();
	for( int j = 0; j < num_iso; j++)
	{

		if(iso_it->second[j].tname != trans_id)continue;

		locate_in_isoform( s, e, iso_it->second[j], vec_match, false);

		vec_all[iso_it->second[j].id].reserve(3);
		vec_all[iso_it->second[j].id].push_back( (uint32_t)flag);
		vec_all[iso_it->second[j].id].push_back(vec_match[6]);	vec_all[iso_it->second[j].id].push_back(vec_match[7]);

		if (max_e < vec_match[0])
		{
			max_e = vec_match[0];
			e_s   = vec_match[2];	e_e   = vec_match[3];
			best_trans = iso_it->second[j].id;
			if( 2 > flag) {flag = 2 + iso_it->second[j].cds_iso;}
		}
		if (max_c < vec_match[1])
		{
			max_c = vec_match[1];
			c_s   = vec_match[4];	c_e   = vec_match[5];
			best_trans = iso_it->second[j].id;
			flag = 4;
		}
	}
		

	if(best_gene == "intergenic")best_trans = best_gene;

	vec_best.clear();
	vec_best.reserve(8);
	vec_best.push_back( (uint32_t)flag);
	vec_best.push_back(max_e);	vec_best.push_back(max_c);
	vec_best.push_back(e_s);	vec_best.push_back(e_e);
	vec_best.push_back(c_s);	vec_best.push_back(c_e);

	return pos;
}
