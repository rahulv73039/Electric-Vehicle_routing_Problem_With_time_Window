//-------------------------------------------------Rahul Verma---(190671)---------------------------------------------------
#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdio>
#include"ilcplex/ilocplex.h"
#include"stdc++.h" 
using namespace std;

struct node {
	double x, y, demand, readyTime, dueDate, serviceTime;  
	string id;
	//char type;
	node(double b, double c, double d ,double e , double f,double g, string h) : x(b) ,y(c), demand(d), readyTime(e), dueDate(f), serviceTime(g) , id(h) {}
	
};

// node-> x;
int main() {
	ifstream fin("input.txt");
	ofstream fout("output.txt");
	string s;
	fin >> s >> s >> s >> s >> s >> s >> s >> s;
	string id;
	map<string, vector<node*>> mp;  
	while (fin >> s) {
		if (s == "Q")break; 
		string a;
		double b, c, d, e, f, g;
		fin >> a >> b >> c >> d >> e >> f >> g;
		node* tmp = new node(b, c, d, e, f, g,s);
		mp[a].push_back(tmp);
	}
	//for (int i = 0; i < mp["c"].size(); i++)cout << mp["c"][i]->x<<endl;
	double Q, C, r, g, v;
	char c;
	fin >> s >> s >> s >> s >> c >> Q >> c >> s >> s >> s >> s >> c >> C >> c >> s >> s >> s >> s >> c >> r >> c >> s >> s >> s >> s >> c >> g >> c >> s >> s >> s >> c >> v >> c;
	//cout << Q << " " << C << " " << r << " " << g << " " << v;
	// All Scanning correct  
	//mp["f'"].push_back(mp["f"][0]); 
	//cout << mp["f"].size() << endl;
	for (int i = 0; i < mp["c"].size(); i++) {
		for (int j = 0; j < mp["f"].size(); j++) {
			mp["f'"].push_back(mp["f"][j]);
		} 
	}  
	
	// by data point of view v'0 = v'n+1 so we declare only one
	mp["v'0"] = mp["d"];
	mp["v'0"].insert(mp["v'0"].end(), mp["c"].begin(), mp["c"].end());
	mp["v'0"].insert(mp["v'0"].end(), mp["f'"].begin(), mp["f'"].end());
	IloEnv env;
	IloModel Model(env);
	int n = mp["v'0"].size();
	IloArray<IloNumVarArray> xij(env,n+1);  
	IloNumVarArray ti(env, n+1, 0, IloInfinity, ILOFLOAT);
	IloNumVarArray yi(env, n + 1, 0, IloInfinity, ILOFLOAT); 
	IloNumVarArray ui(env, n + 1, 0, IloInfinity, ILOFLOAT);
	for (int i = 0; i <n+1 ; i++)xij[i] = IloNumVarArray(env, n+1, 0, 1, ILOINT);
	IloExpr obj(env);  
	vector<vector<double>> dij(n+1,vector<double>(n+1,0)); 
	// calculating dij for  i  in v'0 and  j in v'n+1
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double x1 = mp["v'0"][i]->x, y1 = mp["v'0"][i]->y, x2 = mp["v'0"][j]->x, y2 = mp["v'0"][j]->y; 
			dij[i][j] = sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));

		}
	} 
	// dij for j = n+1
	for (int i = 0; i < n; i++) {
		double x1 = mp["v'0"][i]->x, y1 = mp["v'0"][i]->y, x2 = mp["v'0"][0]->x, y2 = mp["v'0"][0]->y; // since 0 == n+1
	    dij[i][n] = sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));
	}
	for (int j = 0; j < n; j++) {
		double x1 = mp["v'0"][0]->x, y1 = mp["v'0"][0]->y, x2 = mp["v'0"][j]->x, y2 = mp["v'0"][j]->y; // since 0 == n+1
		dij[n][j] = sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));
	}

	for (int i = 0; i < n; i++) {
		for (int j = 1; j < n+1; j++) {
			if(i!=j && !(i==0 && j== n))obj += dij[i][j] * xij[i][j];
		}
	}
	Model.add(IloMinimize(env, obj));
	 // i am assuming 0 is depot 1,mp['c'] size is cust and all till n is f' and n+1 is depot again
	IloRangeArray constr1(env);
 // (2)
	for (int i = 1; i < mp["c"].size() +1; i++) {
		IloExpr exp(env);
		for (int j = 1; j < n+1; j++) {
			if (i != j )exp += xij[i][j];
		}
		constr1.add(exp == 1);
	} 
	//(3)
	for (int i = mp["c"].size() + 1; i < n; i++) {
		IloExpr exp(env);
		for (int j = 1; j < n + 1; j++) {
			if (i != j )exp += xij[i][j];
		}
		constr1.add(exp <= 1);
	} 
	//(4)
	for (int j = 1; j < n; j++) {
		IloExpr exp(env); 
		for (int i = 1; i < n + 1; i++) {
			if (i != j) {
				exp += xij[j][i];
			} 
		}
		for (int i = 0; i < n; i++) {
			if (i != j)exp -= xij[i][j];

		} 
		constr1.add(exp == 0);
	} 
	
	
	//because vehicle is travelling with constant speed say v so tij[i][j] = d[i][j]/v
	vector<vector<double>> tij(n + 1, vector<double>(n + 1, 0));
	for (int i = 0; i < n + 1; i++) {
		for (int j = 0; j < n + 1; j++)tij[i][j] =  dij[i][j]/v;
	} 

	// (5) 
	for (int i = 0; i < mp["c"].size() + 1; i++) {
		for (int j = 1; j < n + 1; j++) {
			// V0 is V'0 from 0 to mp["c"].size+1
			if (i != j) {
				IloExpr exp(env); 
				exp = ti[i] - ti[j] + (tij[i][j] + mp["v'0"][i]->serviceTime) * xij[i][j] - (mp["v'0"][0]->dueDate * (1 - xij[i][j])); 
				constr1.add(exp <= 0);
			}
		}
	}

	//(6)
	for (int i = mp["c"].size() + 1; i < n; i++) {
		for (int j = 1; j < n + 1; j++) {
			if (i != j) {
				IloExpr exp(env); 
				exp = ti[i] - ti[j] + tij[i][j] * xij[i][j] + g * (Q - yi[i]) - (mp["v'0"][0]->dueDate + g * Q) * (1 - xij[i][j]);
				constr1.add(exp <= 0);
			} 
		} 
	} 
	//(7)
	for (int j = 0; j < n; j++) {
		IloExpr exp(env);
		exp = ti[j]; 
		constr1.add(exp >= mp["v'0"][j]->readyTime);
		constr1.add(exp <= mp["v'0"][j]->dueDate);
	}
	// for n+1 
	IloExpr exp1(env);
	exp1 = ti[n];
	constr1.add(exp1 >= mp["v'0"][0]->readyTime);
	constr1.add(exp1 <= mp["v'0"][0]->dueDate);
	//(8) 
	for (int i = 0; i < n; i++) {
		for (int j = 1; j < n + 1; j++) { 
			if (i != j && !(i == 0 && j == n)) {
				IloExpr exp(env);
				exp = ui[i] - ui[j] - mp["v'0"][i]->demand * xij[i][j] + C * (1 - xij[i][j]);
				constr1.add(exp >= 0);
			} 
		}
	} 
	//(9) 
	IloExpr exp(env);
	exp = ui[0];
	constr1.add(exp <= C); 

	//(10)
	// h is charge consumption rate in data r(fuel consumption rate is given) 
	// so h==r taken
	for (int i = 1; i < mp["c"].size() + 1; i++) {
		for (int j = 1; j < n+1; j++) {
			if (i != j) {
				IloExpr exp(env); 
				exp = yi[i] - yi[j] - (r * dij[i][j]) * xij[i][j] + Q * (1 - xij[i][j]);
				constr1.add(exp >= 0);
			}
		}
	} 

	//(11)  
	// here we encouter f'0 we need to add constraint for 0 seperately

	// for i =0 ; and j can not be n here
	for (int j = 1; j < n ; j++) {
			IloExpr exp(env);
			exp = Q - (r * dij[0][j]) * xij[0][j] - yi[j];
			constr1.add(exp >= 0);
	}

	// for i in f'
	for (int i = mp["c"].size() + 1; i < n; i++) {
		for (int j = 1; j < n + 1; j++) {
			if (i != j) {
				IloExpr exp(env); 
				exp = Q - (r * dij[i][j]) * xij[i][j] - yi[j];
				constr1.add(exp >= 0);
			}
		}
	} 
	Model.add(constr1);
	IloCplex mod(Model);
	mod.setOut(env.getNullStream());
	mod.solve(); 
	fout << "Minimum distance" <<" "<<
		mod.getObjValue()<< "\n";
	for (int i = 0; i < n; i++) {
		for (int j = 1; j < n + 1; j++) {
			if (i != j && !(i == 0 && j == n))if (mod.getValue(xij[i][j]) == 1) {
				int tmp;
				if (j == mp["v'0"].size()) tmp = 0;
				else tmp = j;
				if (mp["v'0"][i]->id != mp["v'0"][tmp]->id) {
					
					fout << mp["v'0"][i]->id<<"("<< mod.getValue(ti[i]) << ")" << "---" << tij[i][j] << "---" << mp["v'0"][tmp]->id<< "("<< mod.getValue(ti[j])<<")" << endl;
				}

			}
		}

	}
	
	
	for (int i = 0; i < n; i++) {
		cout << mp["v'0"][i]->id << endl;
	}
	return 0;

}