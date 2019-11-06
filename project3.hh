///////////////////////////////////////////////////////////////////////////////
// maxprotein.hh
//
// Compute the set of foods that maximizes protein, within a calorie budget,
// with the greedy method or exhaustive search.
//
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

#define DEBUG_PERCENT_MATCH 0
#define DEBUG_VAR(var) std::cout << __LINE__ << ":  " << #var << " = " << var << std::endl;

// Simple structure for a single protein
struct Protein
{
	Protein()
	{
		description = "";
		sequence = "";
	}
	Protein(std::string desc, std::string seq)
	{
		description = desc;
		sequence = seq;
	}
	std::string description;
	std::string sequence;
};

// class for BLOSUM penalties.. acts as a matrix holding penalties based
//     on transitions for one amino acid to another
class BlosumPenaltyArray
{
public:
	BlosumPenaltyArray()
	{
		// nothing here
	}
	~BlosumPenaltyArray()
	{
		// nothing here
	}
	BlosumPenaltyArray(BlosumPenaltyArray &that)
	{
		internal_copy(that);
	}
	BlosumPenaltyArray &operator=(BlosumPenaltyArray &that)
	{
		internal_copy(that);
		return *this;
	}

	int get_penalty(char c1, char c2)
	{
		return _penaltyMap[c1][c2];
	}

	void set_penalty(char c1, char c2, int penalty)
	{
		_penaltyMap[c1][c2] = penalty;
	}

	void debug_map()
	{
		for (auto itr1 = _penaltyMap.begin(); itr1 != _penaltyMap.end(); itr1++)
		{
			for (auto itr2 = itr1->second.begin(); itr2 != itr1->second.end(); itr2++)
			{
				std::cout << itr2->second << "  ";
			}
			std::cout << std::endl;
		}
	}

private:
	void internal_copy(BlosumPenaltyArray &that)
	{
		this->_penaltyMap = that._penaltyMap;
	}

	std::map<char, std::map<char, int>> _penaltyMap;
};

// Alias for a vector of shared pointers to Protein objects.
typedef std::vector<std::shared_ptr<Protein>> ProteinVector;

// -------------------------------------------------------------------------
// Load all the proteins from a standard FASTA format file with one line
// per sequence (multi-line sequences are not allowed).
// Returns false on I/O error.
// -------------------------------------------------------------------------
bool load_proteins(ProteinVector &proteins, const std::string &path)
{
	//std::cout << "Loading proteins from [" << path << "]" << std::endl;
	proteins.clear();
	std::ifstream ifs(path.c_str());
	if (!ifs.is_open() || !ifs.good())
	{
		std::cout << "Failed to open [" << path << "]" << std::endl;
		return false;
	}
	int proteinsLoaded = 0;
	bool have_description = false;
	std::shared_ptr<Protein> newProtein = nullptr;
	while (!ifs.eof())
	{
		std::string lineBuffer;
		std::getline(ifs, lineBuffer);
		if (ifs.eof())
		{
			break;
		}
		if (lineBuffer.size() == 0)
		{
			continue;
		}
		if (lineBuffer[0] == '>')
		{
			newProtein = std::shared_ptr<Protein>(new Protein);
			newProtein->description = lineBuffer.substr(1);
			have_description = true;
		}
		else if (have_description)
		{
			newProtein->sequence = lineBuffer;
			proteins.push_back(newProtein);
			proteinsLoaded++;
			have_description = false;
		}
	}

	ifs.close();
	//std::cout << "Loaded " << proteinsLoaded << " proteins from [" << path << "]" << std::endl;

	return true;
}

// -------------------------------------------------------------------------
// Load the BLOSUM penalties from a standard BLOSUM file (matrix format)
// Returns false on I/O error.
// -------------------------------------------------------------------------
bool load_blosum_file(BlosumPenaltyArray &bpa, const std::string &path)
{
	//std::cout << "Loading proteins from [" << path << "]" << std::endl;
	std::ifstream ifs(path.c_str());
	if (!ifs.is_open() || !ifs.good())
	{
		std::cout << "Failed to open [" << path << "]" << std::endl;
		return false;
	}

	std::vector<char> aas; // Create vector to hold our Amino Acids

	while (!ifs.eof())
	{
		std::string lineBuffer;
		std::getline(ifs, lineBuffer);
		if (ifs.eof())
		{
			break;
		}
		if (lineBuffer.size() == 0)
		{
			continue;
		}

		if (lineBuffer[0] == '$')
		{
			std::string buf;
			std::stringstream ss(lineBuffer.substr(1)); // Insert the string into a stream
			while (ss >> buf)
			{
				aas.push_back(buf[0]);
			}
			continue;
		}

		int penalty;
		char thisRowChar = lineBuffer[0];
		std::stringstream ss(lineBuffer.substr(1)); // Insert the string into a stream
		int tokenCount = 0;
		while (ss >> penalty)
		{
			bpa.set_penalty(thisRowChar, aas[tokenCount], penalty);
			tokenCount++;
		}
	}

	//bpa.debug_map();

	return true;
}

// -------------------------------------------------------------------------
int local_alignment(const std::string &string1,
										const std::string &string2,
										BlosumPenaltyArray &bpa,
										std::string &matchString1,
										std::string &matchString2)
{
	int n = string1.size();
    int m = string2.size();
    int D[n+1][m+1];
    char B[n+1][m+1];
    
	//initializing dynamic programming array to zeros
	//and also the backtrack array to '?'
    for (int i = 0 ; i < n + 1 ; i++) {
        for (int j = 0 ; j < m + 1 ; j++) {
            D[i][j] = 0;
            B[i][j] = '?';
        }
    }

    for (int i = 1 ; i < n + 1 ; i++) {
        for (int j = 1 ; j < m + 1; j++) {
            int u = D[i - 1][j] + bpa.get_penalty(string1[i - 1], '*');	//insertion, you are typing something when there shouldn't be anything
            int l = D[i][j - 1] + bpa.get_penalty('*', string2[j - 1]);	//deletion, you didn't type anything when you should've
            int diag = D[i - 1][j - 1] + bpa.get_penalty(string1[i - 1], string2[j - 1]);	//match or substitution
  
			//set the appropriate character to backtracking array regarding the
			//directions l = left, u = up, d = diagonal
            if (l > u) {
                if(l > diag) {
                    B[i][j] = 'l';
                }
                else {
                    B[i][j] = 'd';
                }
            }
            else {
                if (u > diag) {
                    B[i][j] = 'u';
                }
                else {
                    B[i][j] = 'd';
                }
            }

            int temp_max = 0;
			int real_max = 0;

            if (0 > temp_max){
				temp_max = 0;
			}

			temp_max = std::max(u,l);
			real_max = std::max(temp_max, diag);
			D[i][j] = real_max;
        }
    }

    int best_score = 0;
    int bot_i = string1.size();		//last index of string1, meaning it's the bottom row
    int bot_j = 0;					//collumn of which string2 character belongs to

	//get best score
    for (int j = 1 ; j < m + 1 ; j++) {	
        if (D[bot_i][j] > best_score) {
            best_score = D[bot_i][j];
            bot_j = j;
        }
    }

    int i = bot_i;
    int j = bot_j;
    matchString1 = "";
    matchString2 = "";
    
	//Backtracking... start from end and create the new strings 
	//backwards and add '*' if it's a deletion or insertion
	//stops when '?' is reached
	while (B[i][j] != '?'){
        if (B[i][j] == 'u') {							
            matchString1 = matchString1 + string1[i - 1];
            matchString2 = matchString2 + '*';
            i--;
        }
        else if (B[i][j] == 'l') {
            matchString1 = matchString1 + '*';
            matchString2 = matchString2 + string2[j - 1];
            j--;
        }
        else if (B[i][j] == 'd') {
            matchString1 = matchString1 + string1[i - 1];
            matchString2 = matchString2 + string2[j - 1];
            i--;
            j--;
        }
    }

	//After backtracked, backwards strings are finished, reverse the strings
	//to get the proper order that the string was initially passed
	// ex: old backwards string: C->B->A
	// ex: new forwards string: A->B->C
    reverse(matchString1.begin(),matchString1.end());
    reverse(matchString2.begin(),matchString2.end());

	return best_score;
}
// -------------------------------------------------------------------------
std::shared_ptr<Protein> local_alignment_best_match(
		ProteinVector &proteins,
		const std::string &string1,
		BlosumPenaltyArray &bpa,
		std::string &matchString1,
		std::string &matchString2)
{
	std::shared_ptr<Protein> best_protein = proteins[0];
    std::string str1Match = "";
    std::string str2Match = "";
    int best_score = 0;
	int currScore = 0;

	//check each protein's scores and choose the best protein
    for (int i = 0; i < proteins.size(); i++) {
		currScore = local_alignment(string1, (*proteins[i]).sequence, bpa, str1Match, str2Match);
        if (currScore > best_score) {
            best_score = currScore;
            best_protein = proteins[i];
            matchString1 = str1Match;
            matchString2 = str2Match;
        }
    }
	std::cout << "Score:" << best_score << std::endl;
	return best_protein;
}
