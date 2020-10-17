// This is a personal academic project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

// C++ code for 3-site Hamiltonian
#include "itensor/all.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
#include <sstream>
#include <complex>
#include <math.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <exception>
#include "observables.h"
#include "profile.h"
//#include "time_evolution.h"

using namespace itensor;
using namespace std;

//------------------------------------------------------------------
//The function below translate numbers (etc.) into character strings
//the second parameter (optional) is the precision (digits)
template<class T>
inline string to_string(const T &t, unsigned int precision = 0) {
	stringstream ss;
	if (precision > 0)
		ss.precision(precision);
	ss << t;
	return ss.str();
}
//____________________________________________________________
double char2double(char *a) {
	char *end_ptr;
	const double x = strtod(a, &end_ptr);
	if (end_ptr == a || ('\0' != *end_ptr))
		cout << endl << "ERROR :" << a << " is not a valid format for a double."
				<< endl, exit(0);
	return x;
}
//____________________________________________________________
class Parameters: public map<string, double> { // class Parameters have all methods from container "map" (string key, double value) 
public:
	double val(string var_name) const { // .val("C") gives value for parmeter C
		map<string, double>::const_iterator it = find(var_name);
		if (it == end()) {
			cout << "Error: Parameter " << var_name << " is not defined.\n", exit(
					0);
			return 0;
		} else
			return it->second;
	}
	long longval(string var_name) const {
		double v = val(var_name);
		if (abs(double(round(v)) - v) < 1e-6) {
			return long(round(v));
		} else {
			cout << "Error, parameter " << var_name << "=" << v
					<< " is not a long" << endl, exit(0);
			return 0;
		}
	}
	void PRint(ostream &o) const {
		for (map<string, double>::const_iterator it = begin(); it != end();
				it++) {
			o << it->first << "=" << it->second << endl;
		}
	}
	void ReadArguments(int argc, char *argv[]) {
		for (int n = 1; n < argc; n++) {
			string var_name(argv[n]);
			map<string, double>::const_iterator it = find(var_name);

			if (it != end()) {
				n++;
				if (n == argc)
					cerr << "Error: missing value after " << var_name << endl, exit(
							0);
				operator[](var_name) = char2double(argv[n]);
			} else {
				cerr << "Syntax error :" << var_name << endl;
				cout << "List of command-line parameters :";
				PRint (cout);
				exit(0);
			}
		}
	}
};
//_____________________________________________________

// class of PUBLIC parameters
class ThreeSiteParam: public Parameters {
public:
	ThreeSiteParam() { //Constructor
		//Specify below all the allowed parameter names,
		//and their default values
		operator[]("N") = 10; //Length of the chain N-n
		operator[]("begin") = 1;
		operator[]("end") = 10;
		operator[]("J") = 1.0;
		operator[]("tau") = 0.01;  //time step for the unitary evolution
		operator[]("dbeta") = 0.01;
		operator[]("T") = 0;  //Total (final) time
		operator[]("TL") = 100;
		operator[]("TR") = 5;
		operator[]("hL") = 0; //alternating chemical potential or staggered magnetization
		operator[]("hR") = 0;
		operator[]("Entropy") = 0; //entanglement entropy p*log*p between left and right parts of system
		operator[]("Eprof") = 0; // Entropy profile - parameter 0 -> nothing, dt>0 each second=integer parameter
		operator[]("EnergyProf") = 0;
		operator[]("CurrentProf") = 0;
		operator[]("Current") = 0;
		operator[]("Sz") = 0;
		operator[]("SVD_spec") = 0; //SVD spectrum
		operator[]("max_bond") = 4000;  //maximum bond dimension
		operator[]("trunc") = 1e-10;  //maximum truncation error
		operator[]("trunc0") = 1e-10;
		operator[]("energy") = 1e-10;  //convergence criterium on the energy
		operator[]("sweeps") = 999;  //maximum number of sweeps in the DMRG
		operator[]("TrotterOrder") = 2;
		operator[]("antal") = 0;
		operator[]("XXZ") = 0;
		operator[]("PBC") = 0;
		operator[]("beta") = 1;
		operator[]("Energy_beta") = 1;
		operator[]("write_wf") = 0; //0->do not write the w.-f. to disk. if dt>0 => w.-f. to disk (over)written to disk every time=dt*integer.
	}
};
//_____________________________________________________

//I'm creating a 3 site Hamiltonian for system of N sites
class ThreeSiteHamiltonian {
public:
	int dot;
	AutoMPO ampo;  //variable ampo 
	//Define a constructor for this class 'ThreeSiteHamiltonian'
	ThreeSiteHamiltonian(const SiteSet &sites, const ThreeSiteParam &param) :
			ampo(sites) {
		N = length(sites); // size of Hamiltonian comes with onject SiteSet sites, not just as number N
		init(param);   // initializing the Hamiltonian
		dot = N / 2;
		cout << "A Hamiltonian with " << N << " sites was constructed." << endl;
	}
private:
	int N;
	void init(const ThreeSiteParam &param) {    //.init (param)
		const double J = param.val("J");
		const double n = param.val("begin");
		const double hL = param.val("hL");
		const double hR = param.val("hR");
		const double TL = param.val("TL");
		const double TR = param.val("TR");
		//dot = (N) / 2;  //Position of the "dot"
		//cout << "The dot is on site #" << dot << endl;
		//if ((2*N)<=3) cout<<"Error, N="<<N<<" is too small.\n",exit(0);	
		for (int j = n; j <= N - 4; j += 2) {
			// N = number os sites = 2 L
			//Strange coefficients are needed to match with 
			// spin Pauli matrices instead of Sx Sy 
			ampo += J * 4 * 0.25, "S+", j, "S-", j + 4; // 0.5 (SpSm+ SmSp) = SxSx + SySy
			ampo += J * 4 * 0.25, "S-", j, "S+", j + 4;
			ampo += J * -8 * 0.25, "S+", j, "Sz", j + 2, "S-", j + 4;
			ampo += J * -8 * 0.25, "S-", j, "Sz", j + 2, "S+", j + 4; 
		}
		for (int j = n; j <= N; j += 2) {
			double mu;
			if (j <= N / 2) {
				mu = hL * TL;
			} else {
				mu = hR * TR;
			}
			ampo += -J * (mu) * pow(-1, (j + 1) / 2), "Sz", j;
		}
		cout << "H = 3site is construcnted" << endl;
		if (param.val("PBC")) {
			// This part realizes Periodic Boundary Condition (PBC)
			// term (N-1,N,1)
			ampo += J * 4 * 0.25, "S+", N - 3, "S-", 1;
			ampo += J * 4 * 0.25, "S-", N - 3, "S+", 1;
			ampo += J * -8 * 0.25, "S+", N - 3, "Sz", N - 1, "S-", 1;
			ampo += J * -8 * 0.25, "S-", N - 3, "Sz", N - 1, "S+", 1;
			cout << "PBC;\n sites (" << (N - 3 + 1) / 2 << ", "
					<< (N - 1 + 1) / 2 << ", " << (1 + 1) / 2 << "), ";
			// term (N,1,2)
			ampo += J * 4 * 0.25, "S+", N - 1, "S-", 3;
			ampo += J * 4 * 0.25, "S-", N - 1, "S+", 3;
			ampo += J * -8 * 0.25, "S+", N - 1, "Sz", 1, "S-", 3;
			ampo += J * -8 * 0.25, "S-", N - 1, "Sz", 1, "S+", 3;
			cout << "(" << (N - 1 + 2) / 2 << ", " << (1 + 1) / 2 << ", "
					<< (3 + 1) / 2 << ")" << endl;
			cout << "H = 3site is periodic" << endl;
		}
	}
};

//Trotter Gates

class TrotterExp {
public:
	struct TGate {
		int i1 = 0;
		ITensor G;
		TGate() {
		}
		TGate(int i1_, ITensor G_) :
				i1(i1_),  G(G_) {
		}
	};
	TrotterExp(const SiteSet &sites, const ThreeSiteParam &param,
			const complex<double> tau) {
		//N = length(sites);			
		//end = param.val("N");
		initialize(sites, param, tau);
	}
	;
	void initialize(const SiteSet &sites, const ThreeSiteParam &param,
			const complex<double> tau) {
		const int begin = param.val("begin");
		const int end = param.val("end");
		const int order = param.val("TrotterOrder");
		if (order == 1) {
			cout << "trotter 1 scheme" << endl;
			if (imag(tau) < 1e-8) {
				cout << "Temperature evolution, ancillas are not affected"
						<< endl;
				TemperatureGates(begin, end, tau, sites, param);
				TemperatureGates(begin + 2, end, tau, sites, param);
				TemperatureGates(begin + 4, end, tau, sites, param);
			} else {
				TimeGates(begin, end, tau, sites, param);
				TimeGates(begin + 2, end, tau, sites, param);
				TimeGates(begin + 4, end, tau, sites, param);
			}
		} else {
			cout << "trotter 2 scheme" << endl;
			/*
			double a1 = 1. / 6;		// more precise arrpoximation coefficients
			double a2 = 1 - 2. * a1;
			double b1 = (3 - sqrt(3)) / 6.;
			double b2 = 1. / 2 - b1;
			double c1 = 1. / 2;
			*/
			double begin0 = begin; //this variable are needed to change operators ABC
			double begin2 = begin + 2;
			double begin4 = begin + 4;
			//Trotter gates from arxiv.org/abs/1901.04974
			// Eq. (38),(47)
			if (imag(tau) < 1e-8) {
				cout << "Temperature evolution, ancillas are not affected"
						<< endl;
				TemperatureGates(begin0, end, 0.5 * tau, sites, param); //A
				TemperatureGates(begin2, end, 0.5 * tau, sites, param); //B
				TemperatureGates(begin4, end, tau, sites, param); 		//C		
				TemperatureGates(begin2, end, 0.5 * tau, sites, param);  //B
				TemperatureGates(begin0, end, 0.5 * tau, sites, param); //A
				/*
				//this is still a second order but a bit more precise
				TemperatureGates(begin0, end, a1 * tau, sites, param); //A
				TemperatureGates(begin2, end, b1 * tau, sites, param); //B
				TemperatureGates(begin4, end, c1 * tau, sites, param); //C		
				TemperatureGates(begin2, end, b2 * tau, sites, param); //B
				TemperatureGates(begin0, end, a2 * tau, sites, param); //A
				TemperatureGates(begin2, end, b2 * tau, sites, param); //B
				TemperatureGates(begin4, end, c1 * tau, sites, param); //C
				TemperatureGates(begin2, end, b1 * tau, sites, param); //B
				TemperatureGates(begin0, end, a1 * tau, sites, param); //A
				*/
			} else {
				cout << "Time evolutions with Ancilla" << endl;
				TimeGates(begin0, end, 0.5 * tau, sites, param); //A
				TimeGates(begin2, end, 0.5 * tau, sites, param); //B
				TimeGates(begin4, end, tau, sites, param); 	 //C		
				TimeGates(begin2, end, 0.5 * tau, sites, param); //B
				TimeGates(begin0, end, 0.5 * tau, sites, param); //A

				AncillaGates(begin0, end, 0.5 * tau, sites, param); //A
				AncillaGates(begin2, end, 0.5 * tau, sites, param); //B
				AncillaGates(begin4, end, tau, sites, param); 	 //C		
				AncillaGates(begin2, end, 0.5 * tau, sites, param); //B
				AncillaGates(begin0, end, 0.5 * tau, sites, param); //A
				/*
				TimeGates(begin0, end, a1 * tau, sites, param); //A
				TimeGates(begin2, end, b1 * tau, sites, param); //B
				TimeGates(begin4, end, c1 * tau, sites, param); //C		
				TimeGates(begin2, end, b2 * tau, sites, param); //B
				TimeGates(begin0, end, a2 * tau, sites, param); //A
				TimeGates(begin2, end, b2 * tau, sites, param); //B
				TimeGates(begin4, end, c1 * tau, sites, param); //C
				TimeGates(begin2, end, b1 * tau, sites, param); //B
				TimeGates(begin0, end, a1 * tau, sites, param); //A
				*/
			}
		}
	}
	void TemperatureGates(const int begin, const int end,
			const complex<double> tau, const SiteSet &sites,
			const ThreeSiteParam &param) {
		const int step = 6;
		const double J = param.val("J");
		const double hL = param.val("hL");
		const double hR = param.val("hR");
	        const double TL = param.val("TL");
                const double TR = param.val("TR");	
		double mu = 0;
		const int dot = length(sites) / 2;
		cout << "dot in trotter = " << dot << endl;
		//const double hL = param.val("hL");
		//const double hR = param.val("hR");
		//const double TL = param.val("TL");
		//const double TR = param.val("TR");
		for (int j = begin; j < end - 4; j += step) {
			if (j < dot && dot < j + step - 1) {
				cout << "j = [" << j << ", " << j + 2 << ", " << j + 4 << "]"
						<< endl;
				continue; // here we skip part of loop and got to the next j = j+ step
			}
			cout << "j = (" << j << ", " << j + 2 << ", " << j + 4 << ")"
					<< endl;
			auto hh = J * 4 * 0.25
				* op(sites, "Sp", j )
				* op(sites, "Id", j + 2)
				* op(sites, "Sm", j + 4);

			hh += J * 4 * 0.25 
				* op(sites, "Sm", j )
				* op(sites, "Id", j + 2)
				* op(sites, "Sp", j + 4);

			hh += -J * 8 * 0.25
				* op(sites, "Sp", j )
				* op(sites, "Sz", j + 2)
				* op(sites, "Sm", j + 4);

			hh += -J * 8 * 0.25
				* op(sites, "Sm", j ) 
				* op(sites, "Sz", j + 2) 
				* op(sites, "Sp", j + 4);
                        
                        if (j <= dot) {
                        	mu = hL * TL;
                        } else {
                                mu = hR * TR;
                        }
                	hh += -J * mu * pow(-1, (j + 1) / 2) 
				* op(sites, "Sz", j)
				* op(sites, "Id", j + 2)
				* op(sites, "Id", j + 4);
			hh += -J * mu * pow(-1, (j + 1 + 2) / 2)
				* op(sites, "Id", j)
				* op(sites, "Sz", j + 2)
				* op(sites, "Id", j + 4);
			hh += -J * mu * pow(-1, (j + 1 + 4) / 2) 
				* op(sites, "Id", j)
				* op(sites, "Id", j + 2)
				* op(sites, "Sz", j + 4);			
	
			auto G = expHermitian(hh, -tau);
			gates.emplace_back(j, move(G));
		}
	}
	void TimeGates(const int begin, const int end, const complex<double> tau,
			const SiteSet &sites, const ThreeSiteParam &param) {
		const int step = 6;
		const double J = param.val("J");
		for (int j = begin; j < end - 4; j += step) {
			cout << "j = (" << j << ", " << j + 2 << ", " << j + 4 << ")"
					<< endl;
			//this part act on real sites	
			auto hh = J * 4 * 0.25
				* op(sites, "Sp", j )
				* op(sites, "Id", j + 2)
				* op(sites, "Sm", j + 4);

			hh += J * 4 * 0.25 
				* op(sites, "Sm", j )
				* op(sites, "Id", j + 2)
				* op(sites, "Sp", j + 4);

			hh += -J * 8 * 0.25
				* op(sites, "Sp", j )
				* op(sites, "Sz", j + 2)
				* op(sites, "Sm", j + 4);

			hh += -J * 8 * 0.25
				* op(sites, "Sm", j ) 
				* op(sites, "Sz", j + 2) 
				* op(sites, "Sp", j + 4);

			auto G = expHermitian(hh, -tau);
			gates.emplace_back(j, move(G));
		}
	}
	void AncillaGates(const int begin, const int end, const complex<double> tau,
			const SiteSet &sites, const ThreeSiteParam &param) {
		const int step = 6;
		const double J = param.val("J");
		for (int j = begin; j < end - 4; j += step) {
			cout << "j = (" << j << ", " << j + 2 << ", " << j + 4 << ")"
					<< endl;
			//this part act on real sites	
			auto hh = -J * 4 * 0.25
				* op(sites, "Sp", j + 1)
				* op(sites, "Id", j + 3)
				* op(sites, "Sm", j + 5);

			hh += -J * 4 * 0.25 
				* op(sites, "Sm", j + 1)
				* op(sites, "Id", j + 3)
				* op(sites, "Sp", j + 5);

			hh += J * 8 * 0.25
				* op(sites, "Sp", j + 1)
				* op(sites, "Sz", j + 3)
				* op(sites, "Sm", j + 5);

			hh += J * 8 * 0.25
				* op(sites, "Sm", j + 1) 
				* op(sites, "Sz", j + 3) 
				* op(sites, "Sp", j + 5);

			auto G = expHermitian(hh, -tau);
			ancilla_gates.emplace_back(j, move(G));
		}
	}
	void EvolvePhysical(MPS &psi, const Args &args) {
		for (auto &gate : gates) {
			auto j = gate.i1;
			auto &G = gate.G;
			SwapNextSites(psi,j); //swap j,j+1
			SwapNextSites(psi,j+3);//Now physical sites are j+1,j+2,j+3			
			psi.position(j+2);
			auto WF = psi(j + 1) * psi(j + 2) * psi(j + 3);
			WF = G * WF;
			WF /= norm(WF);
			WF.noPrime();
			{
				auto [Uj1,Vj1] = factor(WF, { siteIndex(psi, j + 1), leftLinkIndex(psi, j + 1) }, args);
				auto indR = commonIndex(Uj1, Vj1);
				auto [Uj2,Vj2] = factor(Vj1, { siteIndex(psi, j + 2), indR },args);
				psi.set(j + 1, Uj1);
				psi.set(j + 2, Uj2);
				psi.set(j + 3, Vj2);
				SwapNextSites(psi,j+3);				
				SwapNextSites(psi,j);
			}
		}
	}
	void EvolveAncillas(MPS &psi, const Args &args) {
		for (auto &gate : ancilla_gates) {
			auto j = gate.i1;
			auto &G = gate.G;
			SwapNextSites(psi,j+1); //swap j+1,j+2
			SwapNextSites(psi,j+4);//Now ancilla sites are j+2,j+3,j+4
			psi.position(j+3);
			auto WF = psi(j + 2) * psi(j + 3) * psi(j + 4);
			WF = G * WF;
			WF /= norm(WF);
			WF.noPrime();
			{
				auto [Uj1,Vj1] = factor(WF, { siteIndex(psi, j + 2), leftLinkIndex(psi, j + 2) }, args);
				auto indR = commonIndex(Uj1, Vj1);
				auto [Uj2,Vj2] = factor(Vj1, { siteIndex(psi, j + 3), indR },args);
				psi.set(j + 2, Uj1);
				psi.set(j + 3, Uj2);
				psi.set(j + 4, Vj2);
				SwapNextSites(psi,j+4);				
				SwapNextSites(psi,j+1);
			}
		}
	}
	void Evolve(MPS &psi, const Args &args){
		EvolvePhysical(psi, args);
	//	EvolveAncillas(psi, args);
	}
	void SwapNextSites(MPS &psi, const int j){
		//cout << "in swap" << endl;
		psi.position(j);
		auto WF = psi(j) * psi(j + 1);
		auto [U,V] = factor(WF,
					{ siteIndex(psi, j+1), leftLinkIndex(psi, j) }, {"Truncate=", false});
		psi.set(j, U);
		psi.set(j + 1, V);
		/*		
		auto [U, D, V] = svd(WF,
					{ siteIndex(psi, j+1), leftLinkIndex(psi, j) }, {"Truncate=", false});
		psi.set(j, U);
		psi.set(j + 1, D*V);
		*/
		psi.position(j);
		//psi.orthogonalize();
		//cout << "out swap" << endl;
	}
private:
	vector<TGate> gates;
	vector<TGate> ancilla_gates;
};

void DisconnectChains(MPS &psi, const int j) {
	psi.position(j);
	auto WF = psi(j) * psi(j + 1);
	auto [U, D, V] = svd(WF, { siteIndex(psi, j), leftLinkIndex(psi, j) }, {
			"MaxDim", 1 });
	psi.set(j, U);
	psi.set(j + 1, D * V);
	psi.orthogonalize();
}

// OBSERVABLES are in "observables.h"

// Time evolution e^-iHt = 1-iHt :: REMOVED

//------------------------------------------------------------------
class MyDMRGObserver: public DMRGObserver {
	double previous_energy;
	const double precision;
public:
	MyDMRGObserver(const MPS &psi, double prec = 1e-10) :
			DMRGObserver(psi), precision(prec) {
	}
	bool checkDone(const Args &args = Args::global()) { const
		double energy = args.getReal("Energy", 0);
		cout << "    Energy change:" << energy - previous_energy << endl;
		if (abs(energy - previous_energy) < precision) {
			cout << "   Energy has converged -> stop the DMRG.\n";
			return true;
		} else {
			previous_energy = energy;
			return false;
		}
	}
};

int main(int argc, char *argv[]) {
	LOG_DURATION("MAIN");
	ThreeSiteParam param;
	param.ReadArguments(argc, argv); //Now param contains the parameters, default values or those provided on the command-line

	param.PRint(cout); // Print parameters
	cout.precision(15);
	
	const int N = 2 * param.longval("N");



	SpinHalf sites(N, { "ConserveQNs=", false }); //HILBERT_SPACE = SpinHalf
	MPS psi;
	//--------------------------------------------------------------
	cout << "start defining H" << endl;
	ThreeSiteHamiltonian Ham(sites, param);
	const int dot = Ham.dot;
	auto H = toMPO(Ham.ampo);
	cout << "finish H" << endl;
	;
	auto energy(0);
	psi = MPS(sites);
	cout << "N= " << N << endl;
	for (int i = 1; i < N; i += 2) {// N = 2L, where L is the original system size
									// in this part I'm creating entangled states at sites (j,j+1) 
									// where j is a physical site, j+1 is an ancilla
		auto il = commonIndex(psi(i), psi(i + 1), "Link");//Just to initialize the variable;
		auto ir = commonIndex(psi(i), psi(i + 1), "Link");//Just to initialize the variable;
		auto s1 = sites(i);
		auto s2 = sites(i + 1);
		auto wf = ITensor(); //Just to initialize the variable;
		if (i == 1) {
			ir = commonIndex(psi(i + 1), psi(i + 2), "Link");
			wf = ITensor(s1, s2, ir);
			psi.ref(i) = ITensor(s1);
			psi.ref(i + 1) = ITensor(s2, ir);
			wf.set(s1(1), s2(2), ir(1), ISqrt2); // | +- > - | -+ > state
			wf.set(s1(2), s2(1), ir(1), -ISqrt2);
		} else if (i == N - 1) {
			il = commonIndex(psi(i - 1), psi(i), "Link");
			wf = ITensor(il, s1, s2);
			psi.ref(i) = ITensor(s1, il);
			psi.ref(i + 1) = ITensor(s2);
			wf.set(il(1), s1(1), s2(2), ISqrt2); // | +- > - | -+ > state
			wf.set(il(1), s1(2), s2(1), -ISqrt2);
		} else {
			il = commonIndex(psi(i - 1), psi(i), "Link");
			ir = commonIndex(psi(i + 1), psi(i + 2), "Link");
			wf = ITensor(il, s1, s2, ir);
			psi.ref(i) = ITensor(s1, il);
			psi.ref(i + 1) = ITensor(s2, ir);
			wf.set(il(1), s1(1), s2(2), ir(1), ISqrt2); // | +- > - | -+ > state
			wf.set(il(1), s1(2), s2(1), ir(1), -ISqrt2);
		}
		ITensor D;
		svd(wf, psi.ref(i), D, psi.ref(i + 1)); // put the prepared ancilla "wf" to the |PSI>
		psi.ref(i) *= D;
	}
	//MPS psi2 = psi; 
	energy = inner(psi, H, psi); //<psi|H0|psi>
	cout << "2. Initial energy=" << energy << " .Norm = " << inner(psi, psi)
			<< endl;

	//--------------------------------------------------------------

	double tau = param.val("tau");
	//const int o = param.val("TrotterOrder");
	//MPO expH1, expH2, expH3, expH4, expH5, expH6, expH7;

	auto args = Args("Method=", "DensityMatrix", "Cutoff", param.val("trunc"),
			"MaxDim", param.longval("max_bond"), "Normalize", false);// arguments for time dynamics

	auto args0 = Args("Method=", "DensityMatrix", "Cutoff", param.val("trunc0"),
			"MaxDim", param.longval("max_bond"), "Normalize", false); // arguments for IMAGINARY time == initial temperature state

	// Output . dat files = observables
	ofstream ent, spec, eprof, sz, sz_avrg, spin_cur, spin_cur_prof,
			energy_beta, energy_prof; //here I'm defining output
	ios_base::openmode mode;
	mode = std::ofstream::out; //Erase previous file (if present)

	double dt = param.val("Entropy");
	if (dt != 0) { //Entropy in the center of the chain
		ent.open("Entropy_center.dat", mode);
		ent.precision(15);
		ent << "#time \t Entropy(dot) \t BondDim(dot) \t MaxBondDim\n";
	}
	//---------------------
	dt = param.val("SVD_spec");
	if (dt > 0) { //SVD Spectrum on central bond
		spec.open("SVD_spec.dat", mode);
		spec.precision(15);
		spec << "#Position=" << dot << "\t<SVD_spectrum>\t\ttime\n";
	}
	//---------------------
	dt = param.val("Energy_beta");
	if (dt > 0) { // total Energy vs beta. plus max bond dim at the 3rd column
		energy_beta.open("Energy_beta.dat", mode);
		energy_beta.precision(15);
		energy_beta << "beta \t Energy \n";
	}
	//---------------------
	dt = param.val("Eprof");
	if (dt > 0) { //Full entropy profile
		eprof.open("Entropy_profile.dat", mode);
		eprof.precision(15);
		eprof << "#Position=i-" << dot << setw(16) << "\t Entropy(i)"
				<< setw(16)
				<< "\t Entropy_sqrt \t Entropy_state1 \t time \t\t Bond.Dim(i)\n";
	}
	//---------------------
	dt = param.val("Sz");
	if (dt > 0) { //Full magnetization profile
		sz.open("Sz_profile.dat", mode);
		sz.precision(15);
		sz << "#Position=i-" << "\t<Sz_i>\t" << "\t(-1)^i<Sz_i>\t"
				<< "\t\ttime\n";

		sz_avrg.open("Sz_average_profile.dat", mode);
		sz_avrg.precision(15);
		sz_avrg << "#Position=i-" << "\t0.5<Sz_2+1i> + 0.5<Sz_2+1i>\t"
				<< "\t\ttime\n";

	}
	//---------------------
	dt = param.val("EnergyProf");
	if (dt > 0) { //Energy profile
		energy_prof.open("Energy_profile.dat", mode);
		energy_prof.precision(15);
		energy_prof << "#Position=i-" << "\t<Ham_i>\t" << dot
				<< "\t\ttime(or beta)\n";
	}
	//---------------------
	dt = param.val("CurrentProf");
	if (dt > 0) { //spin current profile
		spin_cur_prof.open("Current_profile.dat", mode);
		spin_cur_prof.precision(15);
		spin_cur_prof << "#Position=i-" << "\t<Cur_i>\t" << dot << "\t\ttime\n";
	}
	//---------------------
	dt = param.val("Current");
	if (dt > 0) { //spin current throw central site
		spin_cur.open("Current.dat", mode);
		spin_cur.precision(15);
		spin_cur << "#time" << "\t<Cur_dot>\t" << "\n";
	}
	//---------------------

	//////////// Trottexp expH //////////
	cout << "Trotter Gates for beta " << endl;
	param["begin"] = 1;
	param["end"] = N;
	const double dbeta = param["dbeta"];
	TrotterExp expH_beta(sites, param, 0.5 * dbeta);
	

	cout << "Trotter Gates Half for beta " << endl;
	param["begin"] = dot + 1;
	TrotterExp expH_beta_half(sites, param, 0.5 * dbeta);

	
        cout << "Trotter Gates for tau" << endl;
        param["begin"] = 1;
	param["hL"] = 0;
        param["hR"] = 0;
	TrotterExp expH(sites, param, Cplx_i * 1.0 * tau);



	//Preparation of biased Temperature state
	const double TL = param.val("TL");
	const double TR = param.val("TR");
	const double beta_min = min(1. / TL, 1. / TR);
	const double beta_max = max(1. / TL, 1. / TR);
	const int beta_steps_min = beta_min / param.val("dbeta");
	const int beta_steps_max = beta_max / param.val("dbeta");
	const double n_steps = param.val("T") / param.val("tau");
	param["begin"] = 1;
	param["hL"] = 0;
	param["hR"] = 0;
	cout << "finish H_half" << endl;
	ThreeSiteHamiltonian Ham_time_dynamic(sites, param);
	auto H_time_dynamic = toMPO(Ham_time_dynamic.ampo);

	const int time_total = beta_steps_max + n_steps;
	for (int n = 0; n <= time_total; ++n) {
		double time = (n - beta_steps_max) * tau; //+param.val("time_shift");

		if (n < beta_steps_max) {
			time = n * dbeta;
			cout << "Beta(1/T) step #" << n << "/" << time_total << "\t beta="
					<< time << endl;
		} else {
			cout << "Time step #" << n << "/" << time_total << "\ttime=" << time
					<< endl;
		}
		cout.flush();
		vector<double> Myspec; //vector which will be the SVD spectrum

		//  ----- Compute observables ---
		//  ----- Entropy between sites i and i+1
		if (param.val("Entropy") != 0) {
			double entr = Entropy(psi, dot, Myspec, 1);
			ent << time << "\t" << setw(16) << setfill('0') << entr << "\t"
					<< BondDim(psi, dot) << "\t" << maxLinkDim(psi) << endl;

			if (param.val("SVD_spec") > 0) {
				if (n % int(param.val("SVD_spec") / tau) == 0) {
					spec << "\"t=" << time << "\"" << endl;
					int si = Myspec.size();
					for (int i = 0; i < si; i++) {
						spec << i + 1 << "\t" << Myspec[i] << "\t" << "\t"
								<< time << endl;
					}
					if (n < time_total)
						spec << endl << endl;
				}
			}
		}
		

		// ------- entanglement Entropy profile -----
		if (param.val("Eprof") > 0)
			if (n % int(param.val("Eprof") / tau) == 0) {
				eprof << "\"t=" << time << "\"" << endl;
				for (int i = 1; i < N; i++) {
					double entr_std = Entropy(psi, i, Myspec, 1); // p log p
					eprof << i + 0.5 - dot << "\t" << setw(16) << setfill('0')
							<< entr_std << "\t" << setw(4) << setfill('0')
							<< BondDim(psi, i) << "\t" << time << endl;
				}
				if (n < time_total)
					eprof << "\n\n"; //I need this part to separate time steps in gnuplot
			}

		// ------- Energy vs 1/Temperature -------
		if (param.val("Energy_beta") > 0) {
			double en = 0;
			double counter = 0;
			int shift = (N / 4) % 2 == 0 ? N / 4 : N / 4 + 1;
			for (int i = dot + 1 - shift; i < dot + 1 + shift; i += 2) {
				en += Energy(psi, sites, i);
				++counter;
			}
			en /= (counter);

			energy_beta << time << "\t"
					<< real(innerC(psi, H, psi)) / real(innerC(psi, psi))
					<< "\t" << en << "\t" << maxLinkDim(psi) << endl;
		}
		// ------- Energy profile -------
		if (param.val("EnergyProf") > 0 && beta_steps_max <= n) {
			if (n % int(param.val("EnergyProf") / tau) == 0) {
				energy_prof << "\"t=" << time << "\"" << endl;
				for (int i = 1; i <= N - 5; i += 2) {
					const double en = Energy(psi, sites, i);
					energy_prof << i / 2 - dot / 2 + 1 << "\t" << en << "\t"
							<< time << endl;
				}

				
				energy_prof << "\n\n"; //I need this part to separate time steps in *.dat files (for gnuplot)
			}
		}

		// ------- Current profile -------
		if (param.val("CurrentProf") > 0) {
			if (n % int(param.val("CurrentProf") / tau) == 0) {
				spin_cur_prof << "\"t=" << time << "\"" << endl;
				for (int i = 1; i <= N - 5; i += 2) {
					const double en = Current(psi, sites, i);
					spin_cur_prof << i / 2 - dot / 2 + 1 << "\t" << en << "\t"
							<< time << endl;
				}

				spin_cur_prof << "\n\n"; //I need this part to separate time steps in *.dat files (for gnuplot)
			}
		}
		// ------- Current throw the central site -------
		if (param.val("Current") > 0) {
			if (n % int(param.val("Current") / tau) == 0) {
				const double cur = Current(psi, sites, dot - 2);
				spin_cur << time << "\t" << cur << endl;
			}
		}

		// ------- Sz profile -------
		if (param.val("Sz") > 0) {
			if (n % int(param.val("Sz") / tau) == 0) {
				sz << "\"t=" << time << "\"" << endl;
				double sz_tot = 0, sz_left = 0, sz_right = 0, sz_dot = 0;
				double sz_odd = 0;
				for (int i = 1; i <= N; i += 2) {
					double s = Sz(psi, sites, i);
					sz_tot += s;
					if (i < dot)
						sz_left += s;
					if (i > dot)
						sz_right += s;
					if (i == dot)
						sz_dot += s;
					sz << i / 2 - dot / 2 + 1 << "\t" << s << "\t"
							<< pow(-1, (i + 1) / 2) * s << "\t" << time << endl;

					if (((i + 1) / 2) % 2 == 1) { //odd site
						sz_odd = s;
					} else {
						sz_avrg << i / 2 - dot / 2 + 1 << "\t"
								<< 0.5 * (s + sz_odd) << "\t" << time << endl;
					}
				}

				{ //I need this part to separate time steps in *.dat files (for gnuplot)
					sz << "\n\n";
					sz_avrg << "\n\n";
				}
				cout << "\n<Sz_left>=" << sz_left << "\t" << "<Sz_right>="
						<< sz_right << "\t" << "<Sz_DOT>=" << sz_dot << "\t"
						<< "<Sz_tot>=" << sz_tot << endl;
			}
		}

		if (n < beta_steps_min) {
			cout << "Temperature evol with H" << endl;
			expH_beta.EvolvePhysical(psi, args0);
			psi.orthogonalize(args);
			//DisconnectChains(psi, dot+1);
			cout << "dot = " << dot + 1 << endl;

			psi.normalize();
		} else if (n < beta_steps_max) {
			cout << "Temperature evol with H_half" << endl;
			expH_beta_half.EvolvePhysical(psi, args0);
			psi.orthogonalize(args);			
			psi.normalize();
		} else {
			cout << "Time evol" << endl;
			expH.Evolve(psi, args);
			psi.orthogonalize(args);
			//psi.normalize();
		}
		cout << "max bond dim = " << maxLinkDim(psi) << endl;
		cout << "Norm = " << real(innerC(psi, psi)) << endl;
		cout << "Energy = " << real(innerC(psi, H_time_dynamic, psi)) << endl
				<< endl;
	}

	cout << "Final observables: \n" << "max bond dim = " << maxLinkDim(psi)
			<< endl << "Energy = " << real(innerC(psi, H_time_dynamic, psi))
			<< endl << endl;

	cout << "\nTime evolution complete.\n Done ! \n";

	return 0;
}

