#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdio>
#include"ilcplex/ilocplex.h"

using namespace std;

int* len_stock = new int(0); // Length of stock
int* n_orders = new int(0); // No. of orders
float* weights; // Order Length
int* demand;   // Demand

ifstream fin("input.txt");
ofstream fout("output.txt");
ofstream lout("logs.txt"); 

/*---------------------------------Initial Pattern Generation-----------------------------------------------------*/
int** initial_pattern() {
	int** pattern = new int* [*n_orders];
	for (int* i = new int(0); *i < *n_orders; (*i)++)
	{
		pattern[*i] = new int[*n_orders];
		for (int* j = new int(0); *j < *n_orders; (*j)++)
		{
			pattern[*i][*j] = 0;
		}
		pattern[*i][*i] = ((*len_stock) / weights[*i]);

	}
	pattern[0][0] = std::floor((*len_stock) / weights[0]);
	lout << "Initial Patterns added to the master problem" << endl;
	for (int* i = new int(0); *i < *n_orders; (*i)++)
	{
		lout << (*i) + 1 << ") [";
		for (int* j = new int(0); *j < *n_orders -1; (*j)++)
		{
			
			 lout << pattern[*i][*j] << " ";
		}
		lout << pattern[*i][(*n_orders) - 1] << "]" << "\n";
	}lout << endl;
	return pattern;
}
/*----------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------LP Formulation--------------------------------------------------------*/
IloModel master_lp_form(IloEnv env, IloNumVarArray x, IloRangeArray constr1 , IloRangeArray constr2, int* iter , int** pattern) {

	IloModel Model(env);

	//Objective function
	IloExpr obj(env);
	for (int* p = new int(0); *p < (*n_orders) + (*iter) - 1; (*p)++) {
		obj += x[*p];
	}
	Model.add(IloMinimize(env, obj));  // Minimisation
	obj.end();


	for (int* i = new int(0); *i < *n_orders; (*i)++)
	{
		IloExpr exp(env);
		for (int* p = new int(0); *p < (*n_orders) + (*iter) - 1; (*p)++)
		{
			exp += (pattern[*p][*i]) * (x[*p]);
		}
		constr1.add(exp >= demand[*i]);
	}

	Model.add(constr1);
	return Model;
}
/*------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------Solving Subproblem(Knapsack)------------------------------------------*/
double Knapsack(double* dual, int* newPattern) {
	IloEnv env;
	IloModel Model(env);
	IloNumVarArray new_pat(env, *n_orders, 0, IloInfinity, ILOINT);

	IloExpr exp0(env);
	for (int* i = new int(0); *i < *n_orders; (*i)++)
	{
		exp0 += dual[*i] * new_pat[*i];
	}
	Model.add(IloMaximize(env, exp0));
	IloExpr exp1(env);
	for (int* i = new int(0); *i < *n_orders; (*i)++)
	{
		exp1 += weights[*i] * new_pat[*i];
	}
	Model.add(exp1 <= (*len_stock));
	IloCplex mod(Model);
	mod.setOut(env.getNullStream());
	if (!mod.solve()) {
		env.error() << "Failed to optimize the Master Problem!!!" << endl;
		throw(-1);
	}
	
	double obj = 0.0;
	obj = mod.getObjValue();
	lout << "Objective value of the sub-problem: " << obj - 1 << endl;
	lout << "New pattern generated : [";
	for (int* i = new int(0); *i < *n_orders; (*i)++)
	{
		newPattern[*i] = mod.getValue(new_pat[*i]);
		 lout << newPattern[*i] << " ";
	}

	lout << newPattern[(*n_orders) - 1] << "]" << endl;
	
	return obj;
}
/*---------------------------------------------------------------------------------------------------*/

int main() {

//------------------Taking Input from file----------------------------------------------------

	char *a = new char[100];
	int* b = new int(0);
	fin>>a>>a>>a>>*b;
	*len_stock = *b;
	fin >> a >> a >> a >> *n_orders;
	fin >> a >> a >> a;
	weights = new(nothrow) float[*n_orders];
	demand = new(nothrow) int[*n_orders];
	for (int* i = new int(0); *i < *n_orders; (*i)++) {
		fin >> weights[*i] >> demand[*i];
	}
//--------------------------------------------------------------------------------------

	int** pattern = initial_pattern();

//------------------------Master Problem-------------------------------------------------
	int** cutting_p = new int* [*n_orders];
	int* npattern = new int[*n_orders];
	int* iter = new int(0);
	double *sub_obj = new double(0.0);
	double *obj_val = new double(0.0);
	double pre_sub_obj = numeric_limits<float>::max();
	while (1) {
		(*iter)++;
		IloEnv env;
		IloNumVarArray x(env, (*n_orders) + (*iter) - 1, 0, IloInfinity, ILOFLOAT);
		IloRangeArray constr1(env);
		IloRangeArray constr2(env);
		IloCplex mod(master_lp_form(env,x, constr1,constr2, iter, pattern));


		mod.setOut(env.getNullStream());
		mod.solve();

		lout << "Iteration " << (*iter) << ": " << endl << "Dual Values: ";
//Dual Value Calculation
		double* dual = new double[*n_orders]();
		for (int* i = new int(0); *i < *n_orders; (*i)++)
		{
			dual[*i] = mod.getDual(constr1[*i]);
			lout << fixed << setprecision(2) << dual[*i] << " ";
		}
		lout << endl;
		*obj_val = mod.getObjValue();
		lout << "Objective value of the master-problem: " << *obj_val << endl;
// generating new pattern from the subproblem (knapsack problem)

		int* newPattern = new int[*n_orders];
		for (int i = 0; i < *n_orders; i++) {
			newPattern[i] = 0;
		}
		*sub_obj = Knapsack(dual, newPattern);


	// check optimality condition or add new pattern
		if ((1 - *sub_obj >= 0) || (pre_sub_obj == *sub_obj)) {
			lout << endl << "Optimality attained for the master-problem.";

			int* index = new int(0);
			for (int i = 0; i < ((*n_orders) + (*iter) - 1); i++) {
				if (ceil(mod.getValue(x[i])) > 0) {
					cutting_p[*index] = pattern[i];
					npattern[*index] = ceil(mod.getValue(x[i]));
					(*index)++;
				}
			}


			break;
		}
		else
		{
			pre_sub_obj = *sub_obj;
			pattern[(*n_orders) + (*iter) - 1] = newPattern;
		}
		lout << endl;

	}

//--------------------------------------------------Output---------------------------------------------------
	fout << "No. of stocks to be cut: " << std::ceil(*obj_val) << endl;
//Waste Calculation----------------------------------------------------------

	float* waste_t = new float(0.0);
	float* waste_i = new float(0.0);
	for (int* p = new int(0); *p < *n_orders; (*p)++) {
		*waste_i = 0;
		for (int* i = new int(0); *i < (*n_orders); (*i)++) {
			*waste_i = (*waste_i) + (weights[*i] * (cutting_p[*p][*i]));
		}
		(*waste_t) += (((*len_stock) - (*waste_i)) * npattern[*p]);
	}
	(*waste_t) = ((*waste_t) / ((*len_stock) * ceil(*obj_val))) ;
	fout << fixed << setprecision(2) << "Waste percentage: " << (*waste_t)*100 << endl;

//------------------------------------------------------------------------------
	fout << "Order Lengths: [";
	for (int* i = new int(0); (*i) < (*n_orders) -1; (*i)++) {
		 fout << fixed << setprecision(2) << weights[*i] << " ";
	}
	fout<< fixed << setprecision(2) << weights[(*n_orders) - 1] << "]" << endl;
//------------------------------------------------------------------------------

	fout << "Cutting Pattern" << "		" << "No.of times cut" << endl;
	for (int* p = new int(0); *p < *n_orders; (*p)++) {
		fout << "[";
		for (int* i = new int(0); *i < (*n_orders) - 1; (*i)++) {
			fout << cutting_p[*p][*i] << " ";
		}
		fout << cutting_p[*p][(*n_orders) - 1] << ']' << "			" << npattern[*p] << endl;
	}
//--------------------------------------------------------------------------------
	return 0;

}