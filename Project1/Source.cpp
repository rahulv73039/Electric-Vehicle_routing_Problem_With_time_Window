//#define _CRT_SECURE_NO_DEPRECATE
#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdio>
#include"ilcplex/ilocplex.h"

using namespace std;

//Declaring variables globally
int* W = new int(0);
int* n_orders = new int(0);
float* weights;
int* d;

ifstream fin("input.txt");
ofstream fout("output.txt");
ofstream lout("logs.txt"); 

int** initial_pattern() {
	int** pattern = new int* [*n_orders];
	for (int* i = new int(0); *i < *n_orders; (*i)++)
	{
		pattern[*i] = new int[*n_orders];
		for (int* j = new int(0); *j < *n_orders; (*j)++)
		{
			pattern[*i][*j] = 0;
		}
		pattern[*i][*i] = ((*W) / weights[*i]);

	}
	pattern[0][0] = std::floor((*W) / weights[0]);
	lout << "Initial Patterns added to the master problem" << endl;
	for (int* i = new int(0); *i < *n_orders; (*i)++)
	{
		lout << (*i) + 1 << ") [";
		for (int* j = new int(0); *j < *n_orders; (*j)++)
		{
			if ((*j) == (*n_orders) - 1)lout << pattern[*i][*j];
			else lout << pattern[*i][*j] << " ";
		}
		lout << "]" << "\n";
	}lout << endl;
	return pattern;
}

IloModel master_lp_form(IloEnv env, IloNumVarArray x, IloRangeArray constr1 , IloRangeArray constr2, int* iter , int** pattern) {

	IloModel Model(env);

	//Objective function: Minimize total waste
	IloExpr obj(env);
	for (int* p = new int(0); *p < (*n_orders) + (*iter) - 1; (*p)++) {
		obj += x[*p];
	}
	Model.add(IloMinimize(env, obj));
	obj.end();


	for (int* i = new int(0); *i < *n_orders; (*i)++)
	{
		IloExpr exp(env);
		for (int* p = new int(0); *p < (*n_orders) + (*iter) - 1; (*p)++)
		{
			exp += (pattern[*p][*i]) * (x[*p]);
		}
		constr1.add(exp >= d[*i]);
	}

	Model.add(constr1);
	return Model;
}

double Knapsack(double* dual, int* newPattern) {
	IloEnv env;
	IloModel Model(env);
	IloNumVarArray newPat(env, *n_orders, 0, IloInfinity, ILOINT);

	IloExpr exp0(env);
	for (int* i = new int(0); *i < *n_orders; (*i)++)
	{
		exp0 += dual[*i] * newPat[*i];
	}
	Model.add(IloMaximize(env, exp0));
	IloExpr exp1(env);
	for (int* i = new int(0); *i < *n_orders; (*i)++)
	{
		exp1 += weights[*i] * newPat[*i];
	}
	Model.add(exp1 <= (*W));
	IloCplex mp(Model);
	mp.setOut(env.getNullStream());
	if (!mp.solve()) {
		env.error() << "Failed to optimize the Master Problem!!!" << endl;
		cout << "fail to optimize the master problem" << endl;
		throw(-1);
	}
	
	double obj = 0.0;
	obj = mp.getObjValue();
	lout << "Objective value of the sub-problem: " << obj - 1 << endl;
	lout << "New pattern generated : [";
	for (int* i = new int(0); *i < *n_orders; (*i)++)
	{
		newPattern[*i] = mp.getValue(newPat[*i]);
		if ((*i) == ((*n_orders) - 1))lout << newPattern[*i];
		else lout << newPattern[*i] << " ";
	}

	lout << "]" << endl;
	
	return obj;
}


int main() {

#pragma region ProblemData

	char *a = new char[100];
	int* b = new int(0);
	fin>>a>>a>>a>>*b;
	*W = *b;
	fin >> a >> a >> a >> *n_orders;
	fin >> a >> a >> a;
	weights = new(nothrow) float[*n_orders];
	d = new(nothrow) int[*n_orders];
	for (int* i = new int(0); *i < *n_orders; (*i)++) {
		fin >> weights[*i] >> d[*i];
	}
#pragma endregion
#pragma region Initial_pattern

	int** pattern = initial_pattern();
	
#pragma endregion
#pragma region Master_Problem
	int** cutting_p = new int* [*n_orders];
	int* npattern = new int[*n_orders];
	int* iter = new int(0);
	double SP_obj;
	double obj_val; 
	double pre_SP_obj = numeric_limits<float>::max();
	while (1) {
		(*iter)++;
		IloEnv env;
		IloNumVarArray x(env, (*n_orders) + (*iter) - 1, 0, IloInfinity, ILOFLOAT);
		IloRangeArray constr1(env);
		IloRangeArray constr2(env);
		IloCplex mp(master_lp_form(env,x, constr1,constr2, iter, pattern));


		mp.setOut(env.getNullStream());
		mp.solve();

		lout << "Iteration " << (*iter) << ": " << endl << "Dual Values: ";
		double* dual = new double[*n_orders]();
		for (int* i = new int(0); *i < *n_orders; (*i)++)
		{
			dual[*i] = mp.getDual(constr1[*i]);
			lout << fixed << setprecision(2) << dual[*i] << " ";
		}
		lout << endl;
		obj_val = mp.getObjValue();
		lout << "Objective value of the master-problem: " << obj_val << endl;
		// generate a new pattern from the subproblem (knapsack problem)

		int* newPattern = new int[*n_orders];
		for (int i = 0; i < *n_orders; i++) {
			newPattern[i] = 0;
		}
		SP_obj = Knapsack(dual, newPattern);


		// check the optimality condition/ add the new pattern
		if ((1 - SP_obj >= 0) || (pre_SP_obj == SP_obj)) {
			lout << endl << "Optimality attained for the master-problem.";

			int* index = new int(0);
			for (int i = 0; i < ((*n_orders) + (*iter) - 1); i++) {
				if (ceil(mp.getValue(x[i])) > 0) {
					cutting_p[*index] = pattern[i];
					npattern[*index] = ceil(mp.getValue(x[i]));
					(*index)++;
				}
			}


			break;
		}
		else
		{
			pre_SP_obj = SP_obj;
			pattern[(*n_orders) + (*iter) - 1] = newPattern;
		}
		lout << endl;

	}


	
#pragma endregion
#pragma region waste
	fout << "No. of stocks to be cut: " << std::ceil(obj_val) << endl;
	float* waste_t = new float(0.0);
	float* waste_i = new float(0.0);
	for (int* p = new int(0); *p < *n_orders; (*p)++) {
		*waste_i = 0;
		for (int* i = new int(0); *i < (*n_orders); (*i)++) {
			*waste_i = (*waste_i) + (weights[*i] * (cutting_p[*p][*i]));
		}
		(*waste_t) += (((*W) - (*waste_i)) * npattern[*p]);
	}
	(*waste_t) = ((*waste_t) / ((*W) * std::ceil(obj_val))) * 100;
	fout << fixed << setprecision(2) << "Waste percentage: " << *waste_t << endl;
#pragma endregion

	fout << "Order Lengths: [";
	for (int* i = new int(0); (*i) < (*n_orders); (*i)++) {
		if ((*i) == ((*n_orders) - 1))fout << fixed << setprecision(2) << weights[*i];
		else fout << fixed << setprecision(2) << weights[*i] << " ";
	}fout << "]" << endl;
	fout << "Cutting Pattern" << "		" << "No.of times cut" << endl;
	for (int* p = new int(0); *p < *n_orders; (*p)++) {
		fout << "[";
		for (int* i = new int(0); *i < (*n_orders) - 1; (*i)++) {
			fout << cutting_p[*p][*i] << " ";
		}
		fout << cutting_p[*p][(*n_orders) - 1] << ']' << "			" << npattern[*p] << endl;
	}
	return 0;

}