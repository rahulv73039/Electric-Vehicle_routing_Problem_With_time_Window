#include<iostream>
#include<fstream>
#include<cmath>
#include"ilcplex/ilocplex.h"

using namespace std;

//Declaring variables globally
int* W = new int(0);
int* I = new int(0);
float* w;
int* d;

ifstream fin("input.txt");
ofstream fout("output.txt");
ofstream logt("logs.txt");

double Knapsack(double* dual, int* newPattern) {
	IloEnv env;
	IloModel Model(env);
	IloNumVarArray newPat(env, *I, 0, IloInfinity, ILOINT);

	IloExpr exp0(env);
	for (int* i = new int(0); *i < *I; (*i)++)
	{
		exp0 += dual[*i] * newPat[*i];
	}
	Model.add(IloMaximize(env, exp0));
	IloExpr exp1(env);
	for (int* i = new int(0); *i < *I; (*i)++)
	{
		exp1 += w[*i] * newPat[*i];
	}
	Model.add(exp1 <= (*W));
	IloCplex mp(Model);
	mp.setOut(env.getNullStream());
	if (!mp.solve()) {
		env.error() << "Failed to optimize the Master Problem!!!" << endl;
		cout << "fail to optimize the master problem" << endl;
		throw(-1);
	}
	logt << "New patttern generated : [";
	for (int* i = new int(0); *i < *I; (*i)++)
	{
		newPattern[*i] = mp.getValue(newPat[*i]);
		if ((*i) == ((*I) - 1))logt << newPattern[*i];
		else logt << newPattern[*i] << " ";
	}

	logt << "]" << endl;
	double obj = 0.0;
	obj = mp.getObjValue();
	return obj;
}

int main() {

#pragma region ProblemData

	char* a = new char('a');
	float* b = new float(0.0);
	int* c = new int(0);
	for (int* i = new int(0); *i < 14; (*i)++)
	{
		fin >> *a;
	}
	fin >> *W;
	for (int* i = new int(0); *i < 12; (*i)++)
	{
		fin >> *a;
	}
	fin >> *I;
	w = new(nothrow) float[*I];
	d = new(nothrow) int[*I];
	for (int* i = new int(0); *i < 17; (*i)++)
	{
		fin >> *a;
	}
	for (int* i = new int(0); *i < *I; (*i)++)
	{
		fin >> *b;
		fin >> *c;
		w[*i] = *b;
		d[*i] = *c;
	}
#pragma endregion
#pragma region Initial_pattern
	int** pattern = new int* [*I];
	for (int* i = new int(0); *i < *I; (*i)++)
	{
		pattern[*i] = new int[*I];
		for (int* j = new int(0); *j < *I; (*j)++)
		{
			pattern[*i][*j] = 0;
		}
		pattern[*i][*i] = ((*W) / w[*i]);

	}
	pattern[0][0] = std::floor((*W) / w[0]);
	logt << "Initial Patterns added to the master problem" << endl;
	for (int* i = new int(0); *i < *I; (*i)++)
	{
		logt << (*i) + 1 << ") [";
		for (int* j = new int(0); *j < *I; (*j)++)
		{
			if ((*j) == (*I) - 1)logt << pattern[*i][*j];
			else logt << pattern[*i][*j] << " ";
		}
		logt << "]" << "\n";
	}logt << endl;

	int** cutting_p = new int* [*I];
	int* npattern = new int[*I];
#pragma endregion
#pragma region Master_Problem
	int* iter = new int(0);
	double SP_obj; double obj_val; double pre_SP_obj = std::numeric_limits<float>::max();;
	while (1) {
		(*iter)++;
		IloEnv env;
		IloModel Model(env);
		IloNumVarArray x(env, (*I) + (*iter) - 1, 0, IloInfinity, ILOFLOAT);
		IloRangeArray constr1(env);
		IloRangeArray constr2(env);

		//Objective function: Minimize total waste
		IloExpr obj(env);
		for (int* p = new int(0); *p < (*I) + (*iter) - 1; (*p)++) {
			obj += x[*p];
		}
		Model.add(IloMinimize(env, obj));
		obj.end();


		for (int* i = new int(0); *i < *I; (*i)++)
		{
			IloExpr exp(env);
			for (int* p = new int(0); *p < (*I) + (*iter) - 1; (*p)++)
			{
				exp += (pattern[*p][*i]) * (x[*p]);              /// changed
			}
			constr1.add(exp >= d[*i]);
		}

		Model.add(constr1);
		IloCplex mp(Model);
		mp.setOut(env.getNullStream());
		mp.solve();

		logt << "Iteration " << (*iter) << ": " << endl << "Dual Values: ";
		double* dual = new double[*I]();
		for (int* i = new int(0); *i < *I; (*i)++)
		{
			dual[*i] = mp.getDual(constr1[*i]);
			logt << fixed << setprecision(2) << dual[*i] << " ";
		}
		logt << endl;
		obj_val = mp.getObjValue();
		logt << "Objective value of the master-problem: " << obj_val << endl;
		// generate a new pattern from the subproblem (knapsack problem)

		int* newPattern = new int[*I];
		for (int i = 0; i < *I; i++) {
			newPattern[i] = 0;
		}
		SP_obj = Knapsack(dual, newPattern);
		if (pre_SP_obj == SP_obj)

			logt << "Objective value of the sub-problem: " << SP_obj - 1 << endl;
		// check the optimality condition/ add the new pattern
		if ((1 - SP_obj >= 0) || (pre_SP_obj == SP_obj)) {
			logt << endl << "Optimality attained for the master-problem.";

			int* index = new int(0);
			for (int i = 0; i < ((*I) + (*iter) - 1); i++) {
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
			pattern[(*I) + (*iter) - 1] = newPattern;
		}
		logt << endl;

	}
	fout << "No. of stocks to be cut: " << std::ceil(obj_val) << endl;
#pragma endregion
#pragma region waste
	float* waste_t = new float(0.0);
	float* waste_i = new float(0.0);
	for (int* p = new int(0); *p < *I; (*p)++) {
		*waste_i = 0;
		for (int* i = new int(0); *i < (*I); (*i)++) {
			*waste_i = (*waste_i) + (w[*i] * (cutting_p[*p][*i]));
		}
		(*waste_t) += (((*W) - (*waste_i)) * npattern[*p]);
	}
	(*waste_t) = ((*waste_t) / ((*W) * std::ceil(obj_val))) * 100;
	fout << fixed << setprecision(2) << "Waste percentage: " << *waste_t << endl;
#pragma endregion

	fout << "Order Lengths: [";
	for (int* i = new int(0); (*i) < (*I); (*i)++) {
		if ((*i) == ((*I) - 1))fout << fixed << setprecision(2) << w[*i];
		else fout << fixed << setprecision(2) << w[*i] << " ";
	}fout << "]" << endl;
	fout << "Cutting Pattern" << "		" << "No.of times cut" << endl;
	for (int* p = new int(0); *p < *I; (*p)++) {
		fout << "[";
		for (int* i = new int(0); *i < (*I) - 1; (*i)++) {
			fout << cutting_p[*p][*i] << " ";
		}
		fout << cutting_p[*p][(*I) - 1] << ']' << "			" << npattern[*p] << endl;
	}
	return 0;

}