// Copyright (c) 2009-2010 Satoshi Nakamoto
// Copyright (c) 2009-2017 The Bitcoin Core developers
// Distributed under the MIT software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#include <pow.h>

#include <arith_uint256.h>
#include <chain.h>
#include <primitives/block.h>
#include <uint256.h>


#include <string>
#include <iostream>
#include <list>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <openga.hpp>
#include <sstream>
#include <time.h>
#include <ctime>
#include <stdio.h>
#include <limits>
#include <random>
#include <map>
#include <algorithm>

using std::string;
using std::cout;
using std::endl;

using namespace std;

fstream myfile, myfile1;
string line, line1;
vector <unsigned long long int> hashrate;
int current_blockInterval;
int current_difficultyInterval;
int global_blockInterval;
int global_difficultyInterval;
int global1_blockInterval;
int global1_difficultyInterval;
int global_counter_height;
long double global_difficulty = 1;
long double global_objective1;
long double global_objective2;
long double global1_objective1;
long double global1_objective2;
long STDEV_OF_MINING_POWER = 100000 ;
long AVERAGE_MINING_POWER = 400000 ;


unsigned int GetNextWorkRequired(const CBlockIndex* pindexLast, const CBlockHeader *pblock, const Consensus::Params& params)
{
    assert(pindexLast != nullptr);
    unsigned int nProofOfWorkLimit = UintToArith256(params.powLimit).GetCompact();

    // Only change once per difficulty adjustment interval
    if ((pindexLast->nHeight+1) % params.DifficultyAdjustmentInterval() != 0)
    {
        if (params.fPowAllowMinDifficultyBlocks)
        {
            // Special difficulty rule for testnet:
            // If the new block's timestamp is more than 2* 10 minutes
            // then allow mining of a min-difficulty block.
            if (pblock->GetBlockTime() > pindexLast->GetBlockTime() + params.nPowTargetSpacing*2)
                return nProofOfWorkLimit;
            else
            {
                // Return the last non-special-min-difficulty-rules-block
                const CBlockIndex* pindex = pindexLast;
                while (pindex->pprev && pindex->nHeight % params.DifficultyAdjustmentInterval() != 0 && pindex->nBits == nProofOfWorkLimit)
                    pindex = pindex->pprev;
                return pindex->nBits;
                
            }
        
        return pindexLast->nBits;
    }

    // Go back by what we want to be 14 days worth of blocks
    int nHeightFirst = pindexLast->nHeight - (params.DifficultyAdjustmentInterval()-1);
    assert(nHeightFirst >= 0);
    const CBlockIndex* pindexFirst = pindexLast->GetAncestor(nHeightFirst);
    assert(pindexFirst);
    runGA();

    return CalculateNextWorkRequired(pindexLast, pindexFirst->GetBlockTime(), params);
}

unsigned int CalculateNextWorkRequired(const CBlockIndex* pindexLast, int64_t nFirstBlockTime, const Consensus::Params& params)
{
    if (params.fPowNoRetargeting)
        return pindexLast->nBits;

    // Limit adjustment step
    int64_t nActualTimespan = pindexLast->GetBlockTime() - nFirstBlockTime;
    if (nActualTimespan < params.nPowTargetTimespan/4)
        nActualTimespan = params.nPowTargetTimespan/4;
    if (nActualTimespan > params.nPowTargetTimespan*4)
        nActualTimespan = params.nPowTargetTimespan*4;

    // Retarget
    const arith_uint256 bnPowLimit = UintToArith256(params.powLimit);
    arith_uint256 bnNew;
    bnNew.SetCompact(pindexLast->nBits);
    bnNew *= nActualTimespan;
    bnNew /= params.nPowTargetTimespan;

    if (bnNew > bnPowLimit)
        bnNew = bnPowLimit;

    return bnNew.GetCompact();
}

bool CheckProofOfWork(uint256 hash, unsigned int nBits, const Consensus::Params& params)
{
    bool fNegative;
    bool fOverflow;
    arith_uint256 bnTarget;

    bnTarget.SetCompact(nBits, &fNegative, &fOverflow);

    // Check range
    if (fNegative || bnTarget == 0 || fOverflow || bnTarget > UintToArith256(params.powLimit))
        return false;

    // Check proof of work matches claimed amount
    if (UintToArith256(hash) > bnTarget)
        return false;

    return true;
}


struct MySolution
{
	int blockInterval;
	int difficultyInterval;

	string to_string() const
	{
		return
			string("{")
			+ "blockInterval:" + std::to_string(blockInterval)
			+ ", difficultyInterval:" + std::to_string(difficultyInterval)
			+ "}";
	}

	
};

struct MyMiddleCost
{
	// This is where the results of simulation
	// is stored but not yet finalized.

	double objective1; //std dev of average block time
	double objective2; //std dev of difficulty
	
};

typedef EA::Genetic<MySolution, MyMiddleCost> GA_Type;
typedef EA::GenerationType<MySolution, MyMiddleCost> Generation_Type;

void init_genes(MySolution& p, const std::function<double(void)> &rnd01)
{
	// rnd01() gives a random number in 0~1
	p.blockInterval = 1 + 599 * rnd01();
	//p.blockInterval = 600;
	p.difficultyInterval = 1 + 4031 * rnd01();
	//p.difficultyInterval = 2016;
}

void init_genes1(MySolution& p, const std::function<double(void)> &rnd01)
{
	// rnd01() gives a random number in 0~1
	//p.blockInterval = 1 + 599 * rnd01();
	p.blockInterval = current_blockInterval;
	p.difficultyInterval = current_difficultyInterval;
	//p.difficultyInterval = 40;
}

int genMiningPower()
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0.0, 1.0);
	double r = dis(gen);
	long minValue = 1;
	long mining_power = (r * STDEV_OF_MINING_POWER + AVERAGE_MINING_POWER);
	return max(mining_power,minValue);
}


bool eval_solution(
	const MySolution& p,
	MyMiddleCost &c)
{
	const int& blockInterval = p.blockInterval;
	const int& difficultyInterval = p.difficultyInterval;

	long double nPowTargetTimespan = blockInterval * difficultyInterval;
	long double new_difficulty = 1;
	long double old_difficulty = global_difficulty;
	//long double old_difficulty = 0.03;
	long double min_difficulty = 1;
	vector <long long int> block_time;
	vector <long double> difficulty_history;
	long long int total_sum = 0;
	long double mean = 0, stdDev = 0;
	

	long double difficulty_sum = 0;
	long double difficulty_mean = 0, difficulty_stdDev = 0;
	int first_block_height = global_counter_height - difficultyInterval + 1 ;

	random_device rd;
	mt19937_64 eng(0327);
	std::normal_distribution<double> distr(1000000000, 30000000); // define the range

	////To calculate individual block time
	for (int i = 0 ; i < (10000); i++)
	{
		double r = distr(eng);
		if (r < 0)
		{
			r = r * (-1);
		}
		long long int individualBlockTime = old_difficulty * (pow(2, 32)) / r; //second
		block_time.push_back(individualBlockTime);
		difficulty_history.push_back(old_difficulty);

		//Sum all block time
		total_sum = total_sum + (individualBlockTime);
		difficulty_sum = difficulty_sum + (old_difficulty);

		//If is time to readjust difficulty
		if (((i+1) % difficultyInterval) == 0)
		{
			unsigned int counter1 = 0;
			unsigned int sum = 0;
			int previousBlockHeight = i + 1 - difficultyInterval;
			for (int a = previousBlockHeight; a < (previousBlockHeight + difficultyInterval); a++)
			{
				sum = sum + block_time[a];
				counter1 = counter1 + 1;
			}

			//limit the adjustment
			if (sum < nPowTargetTimespan / 4)
				sum = nPowTargetTimespan / 4;
			if (sum > nPowTargetTimespan * 4)
				sum = nPowTargetTimespan * 4;
			//calculate new difficulty
			new_difficulty = old_difficulty * (nPowTargetTimespan / sum);

			if (new_difficulty < min_difficulty)
			{
				new_difficulty = min_difficulty;
			}
			old_difficulty = new_difficulty;

		}
	}
	//myfile.open("C:/Users/zihau/Desktop/wtf.csv", ios::out | ios::trunc);
	//for (int i = 0; i < block_time.size(); i++)
	//{
	//	if (myfile.is_open())
	//	{
	//		myfile << (block_time[i]/60) << endl;
	//	}

	//
	//file.close();
	//Calculate standard deviation
	//sum = sum * 100;
	mean = total_sum / block_time.size();
	difficulty_mean = difficulty_sum / difficulty_history.size();

	for (int i = 0; i < block_time.size(); i++)
	{
		stdDev = stdDev + pow(block_time[i] - mean, 2);
	}


	for (int i = 0; i < difficulty_history.size(); i++)
	{
		difficulty_stdDev = difficulty_stdDev + pow(difficulty_history[i] - difficulty_mean, 2);
	}


	c.objective1 = sqrt(stdDev / block_time.size());
	c.objective2 = sqrt(difficulty_stdDev / difficulty_history.size());
	block_time.clear();
	difficulty_history.clear();
	return true; // solution is accepted
}

bool eval_solution1(
	const MySolution& p,
	MyMiddleCost &c)
{
	const int& blockInterval = current_blockInterval;
	const int& difficultyInterval = current_difficultyInterval;

	long double nPowTargetTimespan = blockInterval * difficultyInterval;
	long double new_difficulty = 1;
	long double old_difficulty = global_difficulty;
	//long double old_difficulty = 0.03;
	long double min_difficulty = 1;
	vector <long long int> block_time1;
	vector <long double> difficulty_history1;
	long long int total_sum = 0;
	long double mean = 0, stdDev = 0;


	long double difficulty_sum = 0;
	long double difficulty_mean = 0, difficulty_stdDev = 0;
	int first_block_height = global_counter_height - difficultyInterval + 1;
	random_device rd;
	mt19937_64 eng(0327);
	std::normal_distribution<double> distr(1000000000, 30000000); // define the range
	////To calculate individual block time
	for (int i = 0; i < (10000); i++)
	{
		double r = distr(eng);
		if (r < 0)
		{
			r = r * (-1);
		}
		long long int individualBlockTime = old_difficulty * (pow(2, 32)) / r; //second
		block_time1.push_back(individualBlockTime);
		difficulty_history1.push_back(old_difficulty);

		//Sum all block time
		total_sum = total_sum + (individualBlockTime);
		difficulty_sum = difficulty_sum + (old_difficulty);

		//If is time to readjust difficulty
		if (((i + 1) % difficultyInterval) == 0)
		{
			unsigned int counter1 = 0;
			long long int sum = 0;
			int previousBlockHeight = i + 1 - difficultyInterval;
			for (int a = previousBlockHeight; a < (previousBlockHeight + difficultyInterval); a++)
			{
				sum = sum + block_time1[a];
				counter1 = counter1 + 1;
			}

			//limit the adjustment
			if (sum < nPowTargetTimespan / 4)
				sum = nPowTargetTimespan / 4;
			if (sum > nPowTargetTimespan * 4)
				sum = nPowTargetTimespan * 4;
			//calculate new difficulty
			new_difficulty = old_difficulty * (nPowTargetTimespan / sum);

			if (new_difficulty < min_difficulty)
			{
				new_difficulty = min_difficulty;
			}
			old_difficulty = new_difficulty;

		}
	}
	//myfile.open("C:/Users/zihau/Desktop/wtf.csv", ios::out | ios::trunc);
	//for (int i = 0; i < block_time.size(); i++)
	//{
	//	if (myfile.is_open())
	//	{
	//		myfile << (block_time[i]/60) << endl;
	//	}

	//
	//file.close();
	//Calculate standard deviation
	//sum = sum * 100;
	mean = total_sum / block_time1.size();
	difficulty_mean = difficulty_sum / difficulty_history1.size();

	for (int i = 0; i < block_time1.size(); i++)
	{
		stdDev = stdDev + pow(block_time1[i] - mean, 2);
	}


	for (int i = 0; i < difficulty_history1.size(); i++)
	{
		difficulty_stdDev = difficulty_stdDev + pow(difficulty_history1[i] - difficulty_mean, 2);
	}


	c.objective1 = sqrt(stdDev / block_time1.size());
	c.objective2 = sqrt(difficulty_stdDev / difficulty_history1.size());

	block_time1.clear();
	difficulty_history1.clear();
	return true; // solution is accepted

}

MySolution mutate(
	const MySolution& X_base,
	const std::function<double(void)> &rnd01,
	double shrink_scale)
{
	MySolution X_new;
	bool in_range;
	do {
		in_range = true;
		X_new = X_base;
		X_new.blockInterval += 0.2*(rnd01() - rnd01())*shrink_scale;
		//X_new.blockInterval = 600;
		in_range = in_range && (X_new.blockInterval >= 1 && X_new.blockInterval < 600);
		
		X_new.difficultyInterval += 0.2*(rnd01() - rnd01())*shrink_scale;
		//X_new.difficultyInterval = 2016;
		in_range = in_range && (X_new.difficultyInterval >= 1 && X_new.difficultyInterval < 4032);
		//in_range = in_range && (X_new.difficultyInterval > 39 && X_new.difficultyInterval < 41);
	} while (!in_range);
	return X_new;
}

MySolution crossover(
	const MySolution& X1,
	const MySolution& X2,
	const std::function<double(void)> &rnd01)
{
	MySolution X_new;
	double r;
	r = rnd01();
	X_new.blockInterval = r * X1.blockInterval + (1.0 - r)*X2.blockInterval;
	//X_new.blockInterval = 600;

	r = rnd01();
	X_new.difficultyInterval = r * X1.difficultyInterval + (1.0 - r)*X2.difficultyInterval;
	//X_new.difficultyInterval = 2016;
	return X_new;
}

std::vector<double> calculate_MO_objectives(const GA_Type::thisChromosomeType &X)
{
	return {
		X.middle_costs.objective1,
		X.middle_costs.objective2
	};
}

std::ofstream output_file;

void MO_report_generation(
	int generation_number,
	const EA::GenerationType<MySolution, MyMiddleCost> &last_generation,
	const std::vector<unsigned int>& pareto_front)
{
	(void)last_generation;

	//cout << "Generation [" << generation_number << "], ";
	//cout << "Pareto-Front {";
	//for (unsigned int i = 0; i < pareto_front.size(); i++)
	//{
	//	cout << (i > 0 ? "," : "");
	//	cout << pareto_front[i];
	//}
	//cout << "}" << endl;
}

void save_results(const GA_Type &ga_obj)
{
	std::ofstream output_file;
	output_file.open("C:/Users/zihau/Desktop/SDL2-2.0.9/SDLproject/SDLproject/results.txt");
	output_file << "step" << "\t" << "Objective 1" << "\t" << "Objective 2" << "\t" << "Solution" << "\n";

	std::vector<unsigned int> paretofront_indices = ga_obj.last_generation.fronts[0];
	
	for (unsigned int i : paretofront_indices)
	{
		const auto &X = ga_obj.last_generation.chromosomes[i];
		
		output_file << std::fixed << std::setprecision(30)
			<< i << "\t"
			<< X.middle_costs.objective1 << "\t"
			<< X.middle_costs.objective2 << "\t"
			<< X.genes.to_string() << "\n";

		global_blockInterval = X.genes.blockInterval;
		global_difficultyInterval = X.genes.difficultyInterval;
	}	
	output_file.close();
}

void save_results2(const GA_Type &ga_obj)
{

	std::vector<unsigned int> paretofront_indices = ga_obj.last_generation.fronts[0];

	for (unsigned int i : paretofront_indices)
	{
		double objective1 = 0;
		double objective2 = 0;
		const auto &X = ga_obj.last_generation.chromosomes[i];
		
				
		global_objective1 = X.middle_costs.objective1;
		global_objective2 = X.middle_costs.objective2;
		global_blockInterval = X.genes.blockInterval;
		global_difficultyInterval = X.genes.difficultyInterval;
		
	}
}

//for current parameters
void save_results3(const GA_Type &ga_obj)
{

	std::vector<unsigned int> paretofront_indices = ga_obj.last_generation.fronts[0];

	for (unsigned int i : paretofront_indices)
	{
		double objective1 = 0;
		double objective2 = 0;
		const auto &X = ga_obj.last_generation.chromosomes[i];


		global1_objective1 = X.middle_costs.objective1;
		global1_objective2 = X.middle_costs.objective2;
		global1_blockInterval = X.genes.blockInterval;
		global1_difficultyInterval = X.genes.difficultyInterval;

	}
}


int runGA()
{
	myfile.open("C:/Users/zihau/Desktop/hash-rate-graph-v2.csv", ios::in);
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			hashrate.push_back(stoll(line));
		}
	}
	myfile.close();

	EA::Chronometer timer;
	timer.tic();


	std::ofstream output_file;
	output_file.open("C:/Users/zihau/Desktop/SDL2-2.0.9/SDLproject/SDLproject/intermediate_result.csv",ios::app);
	output_file  << "Objective 1" << "," << "GA_Objective 1" << "," << "Objective 2" << "," << "GA_Objective 2" << "," << "blockInterval" << "," << "difficultyInterval"  << "\n";
	bool testExpression = true;
	genMiningPower();
	while (testExpression)
	{
		int blockInterval = 600;
		int difficultyInterval = 2016;
		int difficultyCounter = 0;
		current_blockInterval = blockInterval;
		current_difficultyInterval = difficultyInterval;
		long double nPowTargetTimespan = blockInterval * difficultyInterval;
		long double new_difficulty = 1;
		long double old_difficulty = 1;
		//long double old_difficulty = 0.03;
		long double min_difficulty = old_difficulty;
		vector <long long int> block_time2;
		vector <double> difficulty_history;
		long long int total_sum = 0;
		long double difficulty_sum = 0;
		
		int fixed_blockInterval = 600;
		int fixed_difficultyInterval = 2016;
		long double fixed_nPowTargetTimespan = fixed_blockInterval * fixed_difficultyInterval;
		long double fixed_newDifficulty = 1;
		long double fixed_oldDifficulty = 1;
		long double fixed_minDifficulty = fixed_oldDifficulty;
		vector <long long int> fixed_blockTime;
		vector <long double> fixed_difficultyHistory;
		long long int fixed_totalSum = 0;
		long double fixed_Mean = 0, fixed_stdDev = 0;
		long double fixed_difficultySum = 0;
		long double fixed_difficultyMean = 0, fixed_difficultystdDev = 0;
		unsigned int fixed_counter1 = 0;
		long long int fixed_Sum = 0;
		int fixed_previousBlockHeight = 0;


		std::random_device rd; // obtain a random number from hardware
		std::mt19937 eng(0327); // seed the generator (1010)
		std::normal_distribution<double> distr(1000000000, 30000000); // define the range
		
		////To calculate individual block time
		for (int i = 0; i < 20160; i++)
		{
			global_counter_height = i;
			
			/*
			std::normal_distribution<double> dist(0.0,1.0);
			double wtf = dist(eng);
			double lower_bound = 0;
			double upper_bound = 1;
			std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
			std::default_random_engine re;
			double a_random_double = unif(re);
			double p = 1.0 / 120000000000;
			for (int i = 0; i < 10; i++)
			{
				double interval = (log(a_random_double) / log(1.0-p))/distr(eng);
			}
			*/

			//set the random number
			double r = distr(eng);
			if (r < 0)
			{
				r = r * (-1);
			}

			//calculate block time (With GA)
			long long int individualBlockTime = old_difficulty * (pow(2, 32)) / r; //second
 			block_time2.push_back(individualBlockTime);
			difficulty_history.push_back(old_difficulty);
			difficultyCounter = difficultyCounter + 1;

			//calculate block time (Without GA)
			long long int fixed_individualBlockTime = fixed_oldDifficulty * (pow(2, 32)) / r;
			fixed_blockTime.push_back(fixed_individualBlockTime);
			fixed_difficultyHistory.push_back(fixed_oldDifficulty);


			//Sum all block time & difficulty (with GA)
			total_sum = total_sum + individualBlockTime;
			difficulty_sum = difficulty_sum + old_difficulty;

			//Sum all block time & difficulty (without GA)
			fixed_totalSum = fixed_totalSum + fixed_individualBlockTime;
			fixed_difficultySum = fixed_difficultySum + fixed_oldDifficulty;

			
			//Without GA
			if (((i + 1) % fixed_difficultyInterval) == 0)
			{
				unsigned int fixed_counter1 = 0;
				long long int fixed_Sum = 0;

				int fixed_previousBlockHeight = i + 1 - fixed_difficultyInterval;
				for (int a = fixed_previousBlockHeight; a < (fixed_previousBlockHeight + fixed_difficultyInterval); a++)
				{
					fixed_Sum = fixed_Sum + fixed_blockTime[a];
					fixed_counter1 = fixed_counter1 + 1;
				}
				
				//limit the adjustment
				if (fixed_Sum < fixed_nPowTargetTimespan / 4)
					fixed_Sum = fixed_nPowTargetTimespan / 4;
				if (fixed_Sum > fixed_nPowTargetTimespan * 4)
					fixed_Sum = fixed_nPowTargetTimespan * 4;
				//calculate new difficulty
				fixed_newDifficulty = fixed_oldDifficulty * (fixed_nPowTargetTimespan / fixed_Sum);

				if (fixed_newDifficulty < fixed_minDifficulty)
				{
					fixed_newDifficulty = fixed_minDifficulty;
				}
				fixed_oldDifficulty = fixed_newDifficulty;

			}


			//If is time to readjust difficulty
			//With GA
			if (((difficultyCounter) % difficultyInterval) == 0)
			{
				
				unsigned int counter1 = 0;
				long long int sum = 0;
				int previousBlockHeight = i + 1 - difficultyInterval;
				for (int a = previousBlockHeight; a < (previousBlockHeight + difficultyCounter); a++)
				{
					sum = sum + block_time2[a];
					counter1 = counter1 + 1;
				}

				//limit the adjustment
				if (sum < nPowTargetTimespan / 4)
					sum = nPowTargetTimespan / 4;
				if (sum > nPowTargetTimespan * 4)
					sum = nPowTargetTimespan * 4;
				//calculate new difficulty
				new_difficulty = old_difficulty * (nPowTargetTimespan / sum);

				if (new_difficulty < min_difficulty)
				{
					new_difficulty = min_difficulty;
				}
				old_difficulty = new_difficulty;
				global_difficulty = old_difficulty;
				
				// To evaluate how the current parameters perform
				GA_Type ga_obj1;
				ga_obj1.problem_mode = EA::GA_MODE::NSGA_III;
				ga_obj1.multi_threading = true;
				ga_obj1.verbose = false;
				ga_obj1.population = 2;
				ga_obj1.generation_max = 1;
				ga_obj1.calculate_MO_objectives = calculate_MO_objectives;
				ga_obj1.init_genes = init_genes1;
				ga_obj1.eval_solution = eval_solution1;
				ga_obj1.mutate = mutate;
				ga_obj1.crossover = crossover;
				ga_obj1.MO_report_generation = MO_report_generation;
				ga_obj1.best_stall_max = 1;
				ga_obj1.elite_count = 1;
				ga_obj1.crossover_fraction = 0.7;
				ga_obj1.mutation_rate = 0.2;
				ga_obj1.solve();
				save_results3(ga_obj1);

				//To evalute the other set of parameters
				GA_Type ga_obj;
				ga_obj.problem_mode = EA::GA_MODE::NSGA_III;
				ga_obj.multi_threading = true;
				ga_obj.verbose = false;
				ga_obj.population = 100; 
				ga_obj.generation_max = 10;
				ga_obj.calculate_MO_objectives = calculate_MO_objectives;
				ga_obj.init_genes = init_genes;
				ga_obj.eval_solution = eval_solution;
				ga_obj.mutate = mutate;
				ga_obj.crossover = crossover;
				ga_obj.MO_report_generation = MO_report_generation;
				ga_obj.best_stall_max = 10;
				ga_obj.elite_count = 10;
				ga_obj.crossover_fraction = 0.7;
				ga_obj.mutation_rate = 0.2;
				ga_obj.solve();
				save_results2(ga_obj);
				
				//If objective from GA is better than the current
				if (global_objective1 < global1_objective1 && global_objective2 < global1_objective2)
				{
					cout << "Current Difficulty interval = " << difficultyInterval << endl;
					cout << "Current Block interval = " << blockInterval << endl;
					blockInterval = global_blockInterval;
					current_blockInterval = blockInterval;
					difficultyInterval = global_difficultyInterval;
					current_difficultyInterval = difficultyInterval;
					cout << "New Difficulty interval = " << difficultyInterval << endl;
					cout << "New block interval = " << blockInterval << endl;
					nPowTargetTimespan = blockInterval * difficultyInterval;
				}
				output_file << std::fixed << std::setprecision(30);
				output_file << global1_objective1 << "," << global_objective1 << "," << global1_objective2 << "," << global_objective2 << "," << blockInterval << "," << difficultyInterval << "\n";
				cout << "Block height = " << global_counter_height << endl;
				cout << "-----------------" << endl;
				difficultyCounter = 0;

			}
			
		}
		//Without GA
		fixed_Mean = fixed_totalSum / fixed_blockTime.size();
		fixed_difficultyMean = fixed_difficultySum / fixed_difficultyHistory.size();
		for (int a = 0; a < fixed_blockTime.size(); a++)
		{
			fixed_stdDev = fixed_stdDev + pow(fixed_blockTime[a] - fixed_Mean, 2);
		}
		for (int a = 0; a < fixed_difficultyHistory.size(); a++)
		{
			fixed_difficultystdDev = fixed_difficultystdDev + pow(fixed_difficultyHistory[a] - fixed_difficultyMean, 2);
		}
		double fixed_objective1 = sqrt(fixed_stdDev / fixed_blockTime.size());
		double fixed_objective2 = sqrt(fixed_difficultystdDev / fixed_difficultyHistory.size());

		//With GA
		long double mean = total_sum / block_time2.size();
		long double difficulty_mean = difficulty_sum / difficulty_history.size();
		long double stdDev = 0;
		long double difficulty_stdDev = 0;
		for (int i = 0; i < block_time2.size(); i++)
		{
			stdDev = stdDev + pow(block_time2[i] - mean, 2);
		}
		for (int i = 0; i < difficulty_history.size(); i++)
		{
			difficulty_stdDev = difficulty_stdDev + pow(difficulty_history[i] - difficulty_mean, 2);
		}
		long double ga_objective1 = sqrt(stdDev / block_time2.size());
		long double ga_objective2 = sqrt(difficulty_stdDev / difficulty_history.size());

		testExpression = false;
		cout << "Block interval = " << blockInterval << endl;
		cout << "Difficulty interval = " << difficultyInterval << endl;
		output_file.close();


		std::ofstream output_file1;
		output_file1.open("C:/Users/zihau/Desktop/SDL2-2.0.9/SDLproject/SDLproject/fixed_results.csv", ios::app);
		//output_file1 << "Objective 1" << "," << "Objective 2" << ","<< "\n";
		output_file1 << std::fixed << std::setprecision(30);
		output_file1 << fixed_objective1 << "," << fixed_objective2 << "\n";
		output_file1.close();

		std::ofstream output_file2;
		output_file2.open("C:/Users/zihau/Desktop/SDL2-2.0.9/SDLproject/SDLproject/block_time.csv");
		output_file2 << "Block number" << "," << "Time" << "," << "\n";
		output_file2 << std::fixed << std::setprecision(30);
		for (int i = 0; i < block_time2.size(); i++)
		{
			output_file2 << i << "," << block_time2[i] << "\n";
		}
		
		std::ofstream output_file3;
		output_file3.open("C:/Users/zihau/Desktop/SDL2-2.0.9/SDLproject/SDLproject/fixed_para_block_time.csv");
		output_file3 << "Block number" << "," << "Time" << "," << "\n";
		output_file3 << std::fixed << std::setprecision(30);
		for (int i = 0; i < block_time2.size(); i++)
		{
			output_file3 << i << "," << fixed_blockTime[i] << "\n";
		}

		std::ofstream output_file4;
		output_file4.open("C:/Users/zihau/Desktop/SDL2-2.0.9/SDLproject/SDLproject/difficulty_history.csv");
		output_file4 << "Block number" << "," << "Difficulty" << "," << "\n";
		output_file4 << std::fixed << std::setprecision(30);
		for (int i = 0; i < difficulty_history.size(); i++)
		{
			output_file4 << i << "," << difficulty_history[i] << "\n";
		}

		std::ofstream output_file5;
		output_file5.open("C:/Users/zihau/Desktop/SDL2-2.0.9/SDLproject/SDLproject/fixed_para_difficulty_history.csv");
		output_file5 << "Block number" << "," << "Difficulty" << "," << "\n";
		output_file5 << std::fixed << std::setprecision(30);
		for (int i = 0; i < fixed_difficultyHistory.size(); i++)
		{
			output_file5 << i << "," << fixed_difficultyHistory[i] << "\n";
		}


		std::ofstream output_file6;
		output_file6.open("C:/Users/zihau/Desktop/SDL2-2.0.9/SDLproject/SDLproject/results.csv", ios::app);
		//output_file6 << "Objective 1" << "," << "Objective 2" << "," << "\n";
		output_file6 << std::fixed << std::setprecision(30);
		output_file6 << ga_objective1 << "," << ga_objective2 << "\n";
		output_file6.close();
	}
	
	
	cout << "The problem is optimized in " << timer.toc() << " seconds." << endl;

	return 0;
}