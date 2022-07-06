#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
//This version can find one mismatch

using namespace std;

int match(string read, vector<string> dict)
{
	//read: NGS read, dict: list of sgRNAs, key: constant seq for finding sgRNA location
	//match length: the length of sgRNA to do alignment
	read.erase(read.find_last_not_of(" \r\n\t")+1);
			
	for (int i=0; i<dict.size();i++) // vector variable \B0\B3\BC\F6  
	{
		if (read==dict[i])
			return i; 
		 
	}
	return -1;
}

int main() {
	vector<string> dictionary_REV; //dictionary > guide rna list, dictionary_REV=reverse compliment guide RNA list
	vector<string> dictionary;	//sgRNA library ref seq
	string entry;
	string read_sequence1;
	string read_sequence2;
	int num_seq=0;

	
	int num_reads=0;
	int perfect_matches_F=0;
	int non_perfect_matches=0;
	int perfect_matches_R=0;
	int non_perfect_matches_R=0;
	int coupling_matches=0;
	int decoupling_matches=0;
	int indel=0;
	int no_indel=0;
	
	ifstream sgRNAs;
	sgRNAs.open("sgRNAs-rev compl.txt");
	
	//Reads barcode list
	if (sgRNAs.is_open())
	{
		while (getline(sgRNAs,entry))
		{
			entry.erase(entry.find_last_not_of(" \r\n\t")+1);	//trims new line
			dictionary_REV.push_back(entry);
			num_seq++;
		}
	}
	sgRNAs.close();
	
		
	sgRNAs.open("sgRNAs.txt");
	
	//Reads barcode list
	if (sgRNAs.is_open())
	{
		while (getline(sgRNAs,entry))
		{
			entry.erase(entry.find_last_not_of(" \r\n\t")+1);	//trims new line
			dictionary.push_back(entry);
		}
	}
	sgRNAs.close();
	
	//declares barcode count matrix
	int counter_matrix[num_seq];
	int indel_matrix[num_seq];
	int no_indel_matrix[num_seq];
	
	for (unsigned i=0; i<num_seq;i++)
	{
		counter_matrix[i]=1;//pseudocount to enable logarithmic transformation of dropouts
		indel_matrix[i]=1;
		no_indel_matrix[i]=1;
	}
	
	int y=0;

	//Reads Fastq file
	ifstream fastq1;
	fastq1.open("T1-1.fastq");
	ifstream fastq2;
	fastq2.open("T1-2.fastq");
	int num=0;
	
	while (getline(fastq1,read_sequence1))
	{
		getline(fastq2,read_sequence2);
		num++;
		if(num%4==2)
		{
			num_reads++;
			int sg1=match(read_sequence2.substr(0,20),dictionary); // reverse read; sgrna
			//int sg2=match(read_sequence1.substr(1,20),dictionary_REV); //forward read; target
			int y=read_sequence1.find("CTCGGTGCC");
			
			
			if(sg1>=0)
			{
				if (y>=36)
				{
					counter_matrix[sg1]++;
					
					if(read_sequence1.substr(y-35,10)==dictionary_REV[sg1].substr(10,10)) // right direction? 46-35=11
					{
						coupling_matches++;
						if(read_sequence1.substr(1,20)==dictionary_REV[sg1])//.substr(0.20))
						{
							no_indel_matrix[sg1]++;
							no_indel++; // how?
						}
				
						else //(read_sequence1.substr(y-48,20)!=dictionary_REV[sg1].substr(0.20)); // nessesary? token? why?
						{
							indel_matrix[sg1]++;
							indel++;
						}
					}
					else
						decoupling_matches++;
				}
			}
			else
				non_perfect_matches++;	
		}
	}
/*
	while (getline(fastq1,read_sequence1))
	{
		getline(fastq2,read_sequence2);
		num++;
		if(num%4==2)
		
		string str1=substr(read_sequence1,-"GCACCGACT",26); // from "GCACCGACT", 26 left characters, is it possible to use charaters as position?
		str2=substr(str1,1,6); // former positions of gRNA (6bp) from file1
		str3=substr(read_sequence2,1,6); // former positions of gRNA (6bp) from file2
		num_reads++;
		
		int sg13=match(str2,dictionary);
		int sg23=match(str3,dictionary);
		if (str2=str3) // no coupling
		
			counter_matrix[sg23+sg13*num_seq]++;
				no_coupling_matches++;
		
		
		
		int sg12=read_sequence1.find("GCACCGACT");
		//int sg22=read_sequene2.find("CTAGTCC")==substr(read_sequence2,55,7));
		num_reads++;
		int sg12=match(read_sequence1,dictionary);
		else if (sg12==40)	//no indel
		{
			
				counter_matrix[sg12*num_seq]++;
				No_indel_matches++;
		}
		else  (sg12!=40)	//indel
		
		 && sg22>=0)
			{
				
				Indel_matches++;
		
			}
	
	fastq1.close();
	fastq2.close();
*/	
	int guides_with_reads=0;
	int guides_no_reads=0;
	ofstream count_table;
	count_table.open("library_count.txt");
	count_table<<"sequence\ttarget count\tindel\tno indel\n";
	for (unsigned int x=0; x<num_seq;x++) //x? how to use it
		count_table<<dictionary[x]<<"\t"<<counter_matrix[x]<<"\t"<<indel_matrix[x]<<"\t"<<no_indel_matrix[x]<<"\n";
	count_table.close();
	

	

	
	ofstream summary;
	summary.open("statistics.txt");
	summary<<"Number of coupling_matches: "<<coupling_matches<<"\n";
	summary<<"Number of decoupling_matches: "<<decoupling_matches<<"\n";
	summary<<"Number of indel: "<<indel<<"\n";
	summary<<"Number of no_indel: "<<no_indel<<"\n";
	summary<<"Number of non perfect matches: "<<non_perfect_matches<<"\n";
	summary<<"Coupling ratio: "<<coupling_matches/(coupling_matches+decoupling_matches)<<"\n";
	summary.close();

	return 0;
}
