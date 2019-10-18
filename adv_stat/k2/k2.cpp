#include <vector>
using namespace std;

/*
int alpha(int* node, int* parents, double *dataset)
{
	



}

int parent_instantiation(vector<int> &parents, vector<bool> &mask, double* dataset, int num_features, int d)
{
	vector<int> instants;
	// parents is passed by reference. It contains the number of values each feature can take on
	num_parents = parents.size();
	num_combinations = 0;
	for (int i=0; i<num_parents; i++)
	{
		num_combinations += int(mask[i])*parents[i];
	}
	for (int j=0; j<num_combinations*num_parents; j++)
	{
		instants[j]=
	}
}


int alpha(int idx, int* parents, double* values, vector<bool> &mask, double* dataset, int num_features, int d)
{
	// inside parents, we have to store the number of values that node can take
	// one can think of a mapping from the actual values to the set [0,num_vals-1]
	// the dataset can then be transformed from the actual data into the set of the indexes that correspond to each particular value


	//calculating the number of parents for node idx
	num_parents = 0;
	for (int i=0; i<mask.size(); i++)
	{
		num_parents += int(mask[i]);
	}

	//calculating the number of possible combinations of parent instatiations
	num_combinations=0;
	for (int i=0; i<num_parents; i++)
	{
		num_combinations += int(mask[i])*parents[i];
	}

	//creating an array that counts how many instantiations there are for each combination
	int combinations[num_combinations];
	int alpha = 1;
	int c;

	for (int i=0; i<parents[idx]; i++)
	{
		for (int j=0; j<num_combinations; j++)
			combinations[j]=0;

		c=1;
		for (int j=0; j<d; j++)
		{
			if (dataset[num_features*j+idx]==i)
			{
				for (int l=0; l<num_parents; l++)
				{
					combinations[]+=1;
				}
			}
		}

	}
}
*/

//========================================================================================================================================//
// THINGS MAKE (HOPEFULLY) SENSE FROM NOW ON...


int factorial(int n)
{
	int c = (n<3) ? 1 : n*factorial(n-1);
	return c;
}

int compute_prod_alpha(int idx, int ri, int* instant, int num_parents, int* parents_idxs, int* dataset, int num_features, int d)
{
	// given an instantiation, does the productory of the alphas calculation. 

	//creating an array that counts how many instantiations there are for each value of the node.
	// it is (ri) slots long, so that for each value k of the node idx one collects how many instantiations are there that are like instant
	int combinations[ri]={0};
	int alpha = 1;
	int val;
	bool check=false;

	for (int i=0; i<d; i++)
	{
		// check is used in order to evaluate if a particular row is instantiated according to instant
		// as soon as just one value does not match, it is set to false and so no update is done on combinations
		check=true;
		val=dataset[i*num_features+idx];

		for (int j=0; j<num_parents; j++)
		{
			if (dataset[i*num_features+parents_idxs[j]]!=instant[j])
			{
				check=false;
				break;
			}
		}

		if (check)
		{
			combinations[val]+=1;
		}
	}
	
	// we directly calculate the whole productory of all the alphas
	for (int i=0; i<ri; i++)
	{
		alpha*=factorial(combinations[i])
	}
	return alpha;
}


int compute_Nij(int idx, int* instant, int num_parents, int* parents_idxs, int* dataset, int num_features, int d)
{
	// given an instantiation, does the Nij calculation. 

	// check as before. However, no need of a combinations array this time.
	bool check=false;
	int N=0;

	for (int i=0; i<d; i++)
	{
		check=true;

		for (int j=0; j<num_parents; j++)
		{
			if (dataset[i*num_features+parents_idxs[j]]!=instant[j])
			{
				check=false;
				break;
			}
		}

		if (check)
			N+=1;
	}
	return N;
}

int* instant(int idx, int num_parents, int* parents, int* parents_idxs, int* dataset, int num_features, int d)
{
/* Given:
- the index of  the node i
- the number of its parents 
- the set of the number of values its parents can take
- the indexes of its parents (equivalent to their names)
- the dataset and everything...
this calculates all the possible flavors of parents instantiations


	creating a vector of all possible instances. It is linear, meaning that a user must take proper bites of it

	computing the dimension of the "cartesian product" matrix. The actual array, however, is dimension*num_parents
	this is because:
	
	Suppose two parents, x1=[0,1] and x2=[0,1,2]. Then:

	|(0,0)	(0,1)	(0,2)|	<= matrix of cartesian products. Not useful really... It has dimension 2x3=6
	|(1,0)	(1,1)	(1,2)|

	|0	0|
	|0	1| <= this is the really useful matrix, which has dimesion 2x3x2=12.
	|0	2|
	|1	0|
	|1	1|
	|1	2|

	*/
	int dimension=1;
	for (int i=0; i<num_parents; i++)
	{
		dimension*=parents[parents_idxs[i]];
	}
	
	//the actual 1D array which stores the "useful matrix" above as a flattened array.
	int instants[dimension*num_parents] = {0};
	
	// this array contains the cumulative product of the number of values of each parent.
	// This is used in order to fill instants in a proper way, mimicking an actual cartesian product way of storing the instantiations.
	int parent_utils[num_parents] = {0};
	for (int i=0; i<num_parents; i++)
	{
		parent_utils[j]=1;
		for (int j=0; j<i; j++)
		{
			parent_utils[i]*=parents[parents_idxs[j]];
		}
	}

	// filling instantiations by rows, with i denoting the rows and j denoting the column (i.e., the feature)
	for (int i=0; i<dimension; i++)
	{
		for (int j=0; j<num_parents; j++)
		{	
			instants[i+j]=(i/parent_utils[j])%(parents[parents_idxs[j]]);
		}	
	}

	return instants;
}

int compute_qi(int num_parents, int* parents, int* parents_idxs)
{
// just a dumb funxtion for calculating the value of qi...

	int dimension=1;
	for (int i=0; i<num_parents; i++)
	{
		dimension*=parents[parents_idxs[i]];
	}
	dimension*=num_parents;
	return dimension;
}

double f(int idx, int num_parents, int* parents_idxs, int* parents, int* dataset, int num_features, int d)
{
/*
Given:
- the index of the node, idx
- the number of its parents
- the indexes of its parents
- the number of values each parents can take
- the dataset and all that mumbo jumbo

this calculates the value of f
*/
	int qi, ri, Nij, alpha;
	int* instants;
	double f=1;
	int single_instance[num_parents] = {0};	

	/*
	calculating:
	-the number of possible instantiations
	-the number of possible values the node idx can take
	-all the possible instances
	*/

	qi = compute_qi(num_parents, parents);
	instances = instant(idx, num_parents, parents, parents_idxs, dataset, num_features, d);
	ri = parents[idx];

	/*
	inside the cycle, at each step:
	-one single instantiation is considered, extracting it from the whole vector
	-Nij, the number of times a particular instance appears
	-the productory of the alphas
	*/


	for (int j=0; j<qi; j++)
	{
		for (int i=0; i<num_parents; i++)
		{
			single_instance[i]=instances[j*num_parents+i]:
		}
		Nij = compute_Nij(idx, single_instance, num_parents, parents_idxs, dataset, num_features, d);
		alpha = compute_prod_alpha(idx, ri, single_instance, num_parents, parents_idxs, dataset, num_features, d);
		// formula 20...
		f*=(factorial(ri-1)/factorial(Nij+ri-1))*alpha;
	}
	return f;
}

double fempty(int idx, int* dataset, int num_features, int d)
{
/*
Calculating the probability when no parents are involved.
Given


*/





}

int size(int* parents, int u)
{
// a stupid function that computes the size of the parents set

	int size = 0;
	for (int i=0; i<u; i++)
	{
		size+=(parents[i]==-1)?0:1;
	}
	return size;
}

void reset_parents(int* parents, int u)
{
// an even stupider function that takes the parent set and sets it to all -1s
	for (int i=0; i<u; i++)
	{
		parents[i]=-1;	
	}
}


void k2alg(int* nodes, int* ordering, int* dataset, int num_features, int d, int u)
{
/* where the actual magic happens
Given:
- a set of nodes (nodes)
- an ordering of the nodes (ordering)
- a dataset, together with its dimensions (dataset, num_features, d)
- the maximum number of allowed parents per node (u)
returns the list of optimal parents for each node
*/
	bool flag; //the equivalent of OKToProceed
	int parents[u]={-1}; //the set of parents. Since it can be at most u elements long, that will be its size.
	double P_old, P_new;
	
	for (int i=0; i<num_features; i++)
	{
		reset_parents(parents, u);
		//adding nodes. From j=1 onwards, we will use f and not fempty
		for (j=1; j<u; j++)
		{


		} 
		
		
		
	}






}






















