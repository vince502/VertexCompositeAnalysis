#include <vector>

struct Particle {
    int id;
    double jd;

    // Constructor for Particle
    Particle(int id, double jd) : id(id), jd(jd) {}
};

struct Hadron {
    int id;
    double jd;
    int id1;
    int id2;
    double jd1;
    double jd2;

    // Constructor for Hadron
    Hadron(int id, double jd, int id1, int id2, double jd1, double jd2)
        : id(id), jd(jd), id1(id1), id2(id2), jd1(jd1), jd2(jd2) {}
};

struct Pair {
    Hadron d1;
    Hadron d2;

    // Constructor for Pair
    Pair(const Hadron& d1, const Hadron& d2) : d1(d1), d2(d2) {}
};

// double dR(double a1, double a2, double b1, double b2){
//     return sqrt( (a1-b1)*(a1-b1) + (a2-b2)*(a2-b2)); 
// };
double dR(double a1, double b1){
    return sqrt( (a1-b1)*(a1-b1)); 
};

void testAlgo(){
    std::vector<Hadron> gens={
        Hadron(421, 1.0, -321, 211, 0.101, 1.197),
        Hadron(-421, 1.5, 321, -211, 0.199, 2.198),
    };

    std::vector<Pair> recos ={
        Pair(Hadron(421, 1.0, -321, 211, 0.1, 1.2), Hadron(-421, 1.5, 321, -211, 0.2, 2.2)),
        Pair(Hadron(-421, 1.0, 321, -211, 0.1, 1.2), Hadron(421, 1.5, -321, 211, 0.2, 2.2)),
    };

    auto pass = [](Hadron a, Hadron b){
        if( dR(a.jd1, b.jd1) < 0.03 ) return true;
        return false;
    };

    for(auto reco : recos){
        for(auto gen : gens){
            if(pass(reco.d1, gen)) cout << "Match" << endl;
            if(pass(reco.d2, gen)) cout << "Match" << endl;
            
        }
    }
    

}