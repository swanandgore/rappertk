#include "SCenergyProvider.h"

#include "misc/verbosity.h"
#include "restraints/EDrestraint.h"

#include "builders/Builder.h"
#include "builders/Builder.h"

#include <algorithm>
using std::find;

#include "geometry/Grid.h"
#include "geometry/functions.h"

#include <set>
using std::set;

SCenergyProvider::SCenergyProvider(vector<Builder*> & chibuilders, vector<vector<float> > & PTS, Grid *gr) {
    chibs = chibuilders;
    pts = &PTS;
    grid = gr;
}

void SCenergyProvider::addMTZinfo(const char *mtzFN, const char *FOlabel, const char *FClabel, const char *PHIlabel, int MAPtype, float ESmin, float ESmax) {
    mtzmap = true;
    mtzfn = mtzFN; folabel = FOlabel; fclabel = FClabel; philabel = PHIlabel;
    maptype = MAPtype;
    esmin = ESmin; esmax = ESmax; esmean = ESmax;
}
void SCenergyProvider::addMAPinfo(const char *mtzFN, float ESmin, float ESmax) {
    mtzmap = false;
    mtzfn = mtzFN;
    esmin = ESmin; esmax = ESmax; esmean = ESmax;
}

// find min of available rotamer energies for this sci, as well as arrange validRotamer[sci] indices in increasing order of energy
void SCenergyProvider::orderOnEself(int sci) {
    float temp; int ti;
    for(int i=0; i < eself[sci].size(); i++)
        for(int k=i+1; k < eself[sci].size(); k++) {
            if(eself[sci][i] <= eself[sci][k]) continue;
                temp = eself[sci][i];
                eself[sci][i] = eself[sci][k];
                eself[sci][k] = temp;
                ti = validRotamers[sci][i];
                validRotamers[sci][i] = validRotamers[sci][k];
                validRotamers[sci][k] = ti;
        }
}

void SCenergyProvider::calcEselfmin(int sci) {
    eselfmin[sci] = eself[sci][0];
    for(int i=0; i < eself[sci].size(); i++)
        if(eselfmin[sci] > eself[sci][i])
            eselfmin[sci] = eself[sci][i];
}

// find score based on self-energies of provided rotamers
float SCenergyProvider::scoreGivenRotamers() {
    float score = 0;
    for(int sci=0; sci < chibs.size(); sci++) {
        EDrestraint edsci(chibs[sci]->getOP(), " ");
        if(mtzmap == true)
            edsci = EDrestraint::makeEDrestraintFromMTZ(chibs[sci]->getOP(), " ", mtzfn.c_str(), folabel.c_str(), fclabel.c_str(), philabel.c_str(), maptype , esmin, esmax, esmean);
        else
            edsci = EDrestraint::makeEDrestraintFromMap(chibs[sci]->getOP(), " ", mtzfn.c_str(), esmin, esmax, esmean);
        score += (-1 * edsci.scoreAround(*pts));
    }
    return score;
}


float SCenergyProvider::selfEn(int sci, int ri) {
    if(eself.size() == 0) { eself.resize(chibs.size()); eselfmin.resize(chibs.size()); }
    if(eself[sci].size() == 0) {
        eself[sci].resize(numRot(sci));
        EDrestraint edsci(chibs[sci]->getOP(), " ");
        if(mtzmap == true)
            edsci = EDrestraint::makeEDrestraintFromMTZ(chibs[sci]->getOP(), " ", mtzfn.c_str(), folabel.c_str(), fclabel.c_str(), philabel.c_str(), maptype , esmin, esmax, esmean);
        else
            edsci = EDrestraint::makeEDrestraintFromMap(chibs[sci]->getOP(), " ", mtzfn.c_str(), esmin, esmax, esmean);
        for(int i=0; i < numRot(sci); i++) {
            ((Builder*)chibs[sci])->buildSample(*pts, validRotamers[sci][i]);
            eself[sci][i] = -1 * edsci.scoreAround(*pts);
            //cout << ((Builder*)chibs[sci])->name() << " "<< sci <<" "<< i <<" "<< eself[sci][i] << endl;
            //eself[sci][i] = ((rand()+0.)/RAND_MAX)*100; // a random score assignment
        }
        calcEselfmin(sci); orderOnEself(sci);
    }
    //cout << "selfen " << sci <<" "<< ri <<" "<< eselfmin[sci] <<" "<< eself[sci][ri] << endl;
    if(sci >= eself.size()) cout << "bad access 4 " << sci <<" "<< eself.size() << endl;
    if(ri < 0) return eselfmin[sci];
    if(ri >= eself[sci].size()) cout << "bad access 5 " << ri <<" "<< eself[sci].size() << endl;
    return eself[sci][ri];
}

int SCenergyProvider::numRot(int sci) {
    if(validRotamers.size() > 0) return validRotamers[sci].size();
    validRotamers.resize(chibs.size());
    for(int i=0; i < chibs.size(); i++) {
        for(int ri=0; ; ri++) {
            if(! ((Builder*)chibs[i])->buildSample(*pts, ri)) break; //{ cout << i <<" breaking " << ri << endl; break; }
            vector<int> & bop = chibs[i]->getOP();
            if(grid->add(bop) == true) {
                for(int bi=0; bi < bop.size(); bi++) grid->remove(bop[bi]);
                validRotamers[i].push_back(ri);
            } //else cout << "rejecting a rotamer due to mainchain clashes" << endl;
        }
    }
    if(validRotamers[sci].size() > 0) selfEn(sci,0); //just to ensure that selfEnergies are calcd now itself
    return validRotamers[sci].size();
}

void SCenergyProvider::calcIntxns(int sci1, int sci2) {
    if(sci1 > sci2) { int t; t = sci1; sci1 = sci2; sci2 = t;} 
    if(ipair.find(sci1) != ipair.end() && ipair[sci1].find(sci2) != ipair[sci1].end())
        return;
    if( epair.find(sci1) == epair.end() )
        epair[sci1] = map<int, map<int, map<int, float> > >();
    if( epair[sci1].find(sci2) == epair[sci1].end() )
        epair[sci1][sci2] = map<int, map<int, float> >();
    else return;
    //    if(verbose(6)) cout << "calcIntxn " << sci1 <<" "<< sci2 << endl;
    int nr1 = numRot(sci1);
    int nr2 = numRot(sci2);
    vector<int> & bop1 = chibs[sci1]->getOP();
    vector<int> & bop2 = chibs[sci2]->getOP();
    for(int i=0; i < nr1; i++) {
        int ri1 = validRotamers[sci1][i];
        ((Builder*)chibs[sci1])->buildSample(*pts, ri1);
        if(grid->add(bop1) == false) { cout << "incorrect" << endl; exit(0); }
        for(int k=0; k < nr2; k++) {
            int ri2 = validRotamers[sci2][k];
            ((Builder*)chibs[sci2])->buildSample(*pts, ri2);
            if(epair[sci1][sci2].find(i) == epair[sci1][sci2].end())
                epair[sci1][sci2][i] = map<int,float>();
            epair[sci1][sci2][i][k] = 0;
            if( ! grid->add(bop2) ) {
                epair[sci1][sci2][i][k] = 1e10;
		//                if(verbose(6)) cout << "CLASH " << sci1<<" "<<sci2<<" "<< i <<" "<< k << " " << epair[sci1][sci2][i][k] << endl;
            } else
                for(int bi2=0; bi2 < bop2.size(); bi2++) grid->remove(bop2[bi2]);
        }
        for(int bi1=0; bi1 < bop1.size(); bi1++) grid->remove(bop1[bi1]);
    }
    calcIEpairmin(sci1, sci2);
}

void SCenergyProvider::calcIEpairmin(int sci1, int sci2) {
    if(sci1 > sci2) { int t; t = sci1; sci1 = sci2; sci2 = t; }
    if( ipair.find(sci1) == ipair.end() ) {
        ipair[sci1] = map<int,bool>();
        epairmin[sci1] = map<int,float>();
    }
    if( ipair[sci1].find(sci2) == ipair[sci1].end() ) {
        ipair[sci1][sci2] = false;
        epairmin[sci1][sci2] = 0;
    }
    epairmin[sci1][sci2] = epair[sci1][sci2][0][0];
    ipair[sci1][sci2] = false;
    for(int i=0; i < numRot(sci1); i++) {
        for(int k=0; k < numRot(sci2); k++) {
            if(epair[sci1][sci2][i][k] > 0) ipair[sci1][sci2] = true;
            if(epairmin[sci1][sci2] > epair[sci1][sci2][i][k])
                epairmin[sci1][sci2] = epair[sci1][sci2][i][k];
        }
    }
    if(ipair[sci1][sci2] == false) epair[sci1].erase(sci2);
}


float SCenergyProvider::pairEn(int sci1, int sci2, int ri1, int ri2) {
    if(sci1 > sci2) {
        int t; t = sci1; sci1 = sci2; sci2 = t;
        t = ri1; ri1 = ri2; ri2 = t;
    }
    calcIntxns(sci1,sci2);
    if(! interact(sci1,sci2)) return 0;
    if( epair.find(sci1) == epair.end() || epair[sci1].find(sci2) == epair[sci1].end() || sci1 < 0 || sci2 < 0 )
        cout << "bad access 1 " << sci1<<" "<<sci2<<" "<<ri1<<" "<<ri2<<" "<< epair.size() <<" "<< epair[sci1].size() << endl;
    if( (ri1 < 0 && ri2 >= 0) || (ri2 < 0 && ri1 >= 0) )
        cout << "bad access 3 " << sci1<<" "<<sci2<<" "<<ri1<<" "<<ri2<<" "<< epair.size() <<" "<< epair[sci1].size() << endl;
    if(ri1 < 0 && ri2 < 0) return epairmin[sci1][sci2];
    if(ri1<0 || ri2<0 || ri1 >= epair[sci1][sci2].size() || ri2 >= epair[sci1][sci2][ri1].size())
        cout << "bad access 2 " << sci1<<" "<<sci2<<" "<<ri1<<" "<<ri2<<" "<< epair.size() <<" "<< epair[sci1].size() << endl;
    return epair[sci1][sci2][ri1][ri2];
}

bool SCenergyProvider::interact(int sci1, int sci2) {
    if(sci1 > sci2) { int t; t = sci1; sci1 = sci2; sci2 = t;} 
    calcIntxns(sci1,sci2);
    return ipair[sci1][sci2];
}

// assume that rotamer states have been assigned from SCIs[0...j] both including
// estimate the best possible energy obtainable from j+1...N within this context
float SCenergyProvider::findEbound(vector<int> & SCIs, int J) {
    float E = 0;
    for(int i=J+1; i < SCIs.size(); i++) E += selfEn(SCIs[i]);
    for(int i=J+1; i < SCIs.size(); i++)
        for(int k=0; k <= J; k++) E += pairEn(SCIs[i],SCIs[k]);
    return E;
}


// assume that rotamer states have been assigned from SCIs[0...j] both including
// find the SIGMA Eself + SIGMA Epair for that assignment
float SCenergyProvider::findEsofar(vector<int> & SCIs, int J, vector<int> & assign, bool verbose) {
    for(int i=0; i <= J; i++)
        if(assign[i] < 0 || assign[i] >=  validRotamers[SCIs[i]].size()) {
            cout << "oops111 " << SCIs[i] <<" "<< assign[i]<< endl; exit(0);
        }
    float E = 0;
    for(int i=0; i <= J; i++) {
        E += selfEn(SCIs[i],assign[i]);
	//        if(verbose) cout << "seeE " << SCIs[i] <<" "<< assign[i] <<" "<< selfEn(SCIs[i], assign[i]) <<" "<< validRotamers[SCIs[i]][assign[i]] <<" "<< numRot(SCIs[i]) << endl;
    }
    for(int i=0; i <= J; i++)
        for(int j=i+1; j <= J; j++) {
        //    cout << "findEsofar " << SCIs[i] <<" "<< SCIs[j] <<" "<< i <<" "<< j <<" "<< SCIs.size() << endl;
            E+= pairEn( SCIs[i], SCIs[j], assign[i], assign[j] ) ;
	    //            if(verbose) cout << "seeE " << SCIs[i]<<" "<<SCIs[j]<<" "<< assign[i]<<" "<< assign[j]<<" "<< pairEn( SCIs[i], SCIs[j], assign[i], assign[j] ) << endl;
        }
    return E;
}

// SCIs is a components containing its articulation point A
// for each rotameric state of A, find best possible assignments to SCIs
// modify A's self energy and artiAssig accordingly
void SCenergyProvider::collapse(vector<int> SCIs, int A) {
    Ebest = 0;
    for(int i=0; i < SCIs.size(); i++) // sort acc to num-rot, but keep A as 0th in SCIs
        for(int j=i+1; j < SCIs.size(); j++)
            if(numRot(SCIs[i]) > numRot(SCIs[j])) {
                int t = SCIs[i];
                SCIs[i] = SCIs[j];
                SCIs[j] = t;
            }
    int numposs = 1;
    for(int i=0; i < SCIs.size(); i++) numposs *= numRot(SCIs[i]);
    int iA = 0;
    for(; iA < SCIs.size(); iA++) if(SCIs[iA]==A) break;
    for(int i=iA-1; i >= 0; i--) SCIs[i+1] = SCIs[i];
    SCIs[0] = A;
    cout << "COLLAPSING component " << SCIs.size() <<" (";
    for(int i=0; i < SCIs.size(); i++) cout << SCIs[i] <<" ";
    cout << ") --- into " << A << " -------------------" << endl;

    artiOrder.push_back(A);
    artiAssig.push_back( vector< map<int,int> > (numRot(A)) );
    artiEn.push_back( vector<float> (numRot(A)) );
    for(int ri=0; ri < numRot(A); ri++) {
        cout << "ARTI " << A <<" ^^^^^ rot "<< ri <<" of "<< numRot(A) << endl;
        bestAssig = vector<int>(SCIs.size(),0);
        bestAssig[0] = ri;
        vector<int> currAssig = bestAssig;
        Ebest = findEsofar(SCIs, SCIs.size()-1, bestAssig);
	/*        if(verbose(6)) {
            cout << "Estart " << Ebest <<"--  "<< numposs <<"  (";
            for(int i=0; i < bestAssig.size(); i++) cout << bestAssig[i] <<" "; cout <<") (";
            for(int i=0; i < SCIs.size(); i++) cout << numRot(SCIs[i]) <<" "; cout <<")"<< endl;
        }
	*/

        assign(SCIs, 1, currAssig);
        // now this updates the self-energy of arti-pt
        eself[A][ri] = Ebest;
        artiEn[artiEn.size()-1] [ri] = Ebest;
        for(int i=0; i < SCIs.size(); i++)
            artiAssig[artiAssig.size()-1] [ri] [SCIs[i]] = bestAssig[i];
    }
    // update eselfmin too
    calcEselfmin(A);
    cout << "collapse finished " << eselfmin[A] << endl;
    if(eselfmin[A] > 1e9) {
        cout << "componnet cudnt be solved" << endl;
        exit(0);
    }
}

//assign(j)
//    for rotamer r of sidechain j,
//        assign r to j
//        if Esofar(0..j) is bad, continue
//        if j == last and Esofar < Ebest, replace 
//        find Ebound(0..j, j+1..)
//        if j is last
//            if Esofar < Ebest, update best-assignment-so-far
//        else
//            if Ebest < Esofar + Ebound, continue
//            else assign(j+1)
void SCenergyProvider::assign(vector<int> & SCIs, int J, vector<int> & currAssig) {
    //cout << "ASSIGN " << J << endl;
    for(int ri=0; ri < numRot(SCIs[J]); ri++) {
        currAssig[J] = ri;
        float Esofar = findEsofar(SCIs, J, currAssig);
        if(J == SCIs.size()-1) {
	  //if(verbose(6)) cout <<"Ecompare-last " << Esofar <<" "<< Ebest <<" (";
	  // if(verbose(6)) { for(int i=0; i < currAssig.size(); i++) cout << currAssig[i] <<" "; cout <<")"<< endl; }
            if(Esofar < Ebest) { /*cout <<"Echange " << Esofar <<" "<< Ebest << endl;*/ Ebest = Esofar; bestAssig = currAssig; }
            continue;
        }
        if(Esofar > 1e9) continue;
        float Ebound = findEbound(SCIs, J);
        if(Ebest < Ebound + Esofar) {
	  //            if(verbose(6)) cout << "Eshort! " <<J<<" "<< Ebest <<" "<< Ebound <<" "<< Esofar << endl;
            continue;
        }
	//        else { if(verbose(6)) cout << "Enext " <<J<<" "<< Ebest <<" "<< Ebound <<" "<< Esofar << endl; }
        if(J+1 < SCIs.size()) assign(SCIs, J+1, currAssig);
    }
}

// after all components are collapsed in some order (as in bestAssig),
// global minimal assignment is to be found by travering in reverse direction
void SCenergyProvider::join(vector<int> & finalAssig) {
    // construct components collapsed so far using artiAssig
    finalAssig = vector<int> (chibs.size(), -1);
    vector<vector<int> > components (artiAssig.size()), comps;
    vector<bool> ignored( chibs.size(), true );
    for(int i=0; i < artiAssig.size(); i++)
        for(map<int,int>::iterator it = artiAssig[i][0].begin(); it != artiAssig[i][0].end(); ++it) {
            components[i].push_back(it->first);
            ignored[it->first] = false;
        }
    int nign = 0;
    for(int i=0; i < ignored.size(); i++) {
        if(! ignored[i]) continue;
        nign++;
        int ri = 0;
        for(int k=0; k < numRot(i); k++)
            if(selfEn(i,ri) > selfEn(i,k)) ri = k;
        finalAssig[i] = ri;
        cout << "node-alone " << i <<" chose "<< ri <<" of "<< numRot(i) << " energy " << selfEn(i,ri) << endl;
    }
    cout << "Num-nodes-alone " << nign << endl;
    // go thru components in reverse order. find when anchoring is necessary
    for(int ci=components.size()-1; ci >= 0; ci --) {
        int A = artiOrder[ci];
        bool findAnchor = true;
        //if(ci < components.size()-1)
        //    for(int i=0; i < components[ci+1].size(); i++)
        //        if( find(components[ci].begin(), components[ci].end(), components[ci+1][i]) != components[ci].end() )
        //            { findAnchor = false; break; }
        if(finalAssig[A] < 0) findAnchor = true;
        else findAnchor = false;

        if(findAnchor) cout << "findAnchor" << endl;
        cout << "("; for(int i=0; i < components[ci].size(); i++) cout << components[ci][i] <<" ";
        cout << ") into " << A << endl;

        int ri = 0;
        if(findAnchor) {// choose the best energy rotamer if reqd
            //for(int i=0; i < numRot(A); i++) if(selfEn(A,ri) > selfEn(A,i)) ri = i;
            for(int i=0; i < artiEn[ci].size(); i++) if(artiEn[ci][ri] > artiEn[ci][i]) ri = i;
        }
        else ri = finalAssig[A]; // or use already assigned rotamer
        if(ri < 0) { cout << "ri error " << A<<" "<<ri << endl; exit(0); }
        finalAssig[A] = ri;
        for(map<int,int>::iterator it = artiAssig[ci][ri].begin(); it != artiAssig[ci][ri].end(); ++it)
            if(it->first != A) finalAssig[it->first] = it->second;
    }
}

void SCenergyProvider::resetEnergy() {
    validRotamers.clear();
    eself.clear();
    eselfmin.clear();

    epair.clear();
    epairmin.clear();
    ipair.clear();

    artiOrder.clear();
    artiAssig.clear();
    artiEn.clear();

    bestAssig.clear();
}

void SCenergyProvider::buildAssig(vector<int> & assig) {
    for(int i=0; i < chibs.size(); i++) {
        int ri = validRotamers [i] [ assig[i] ];
        ((Builder*)chibs[i])->buildSample(*pts, ri);
    }
}

// make random assignment and calc energy for it
void SCenergyProvider::energyRandAssig() {
    vector<int> SCIs, randAssig;
    for(int i=0; i < validRotamers.size(); i++) {
        randAssig.push_back(int( floor((rand()/(RAND_MAX+0.)) * numRot(i)) ));
        SCIs.push_back(i);
    }
    float Esofar = findEsofar(SCIs, SCIs.size()-1, randAssig);
    cout << "Erand " << Esofar << endl;
}

void SCenergyProvider::exhaustiveEnergy() {
    vector<int> assig(chibs.size());
    exEn(assig, 0);
}

void SCenergyProvider::exEn(vector<int> & assig, int pos) {
    if(pos == chibs.size()) {
        vector<int> SCIs;
        for(int k=0; k < assig.size(); k++) SCIs.push_back(k);
        float E = findEsofar(SCIs, SCIs.size()-1, assig);
        cout << "Eexha " << E <<" ";
        for(int k=0; k < assig.size(); k++) cout << assig[k]<<":"<< validRotamers[SCIs[k]][assig[k]]<<" ";
        cout << endl;
        return;
    }
    for(int i=0; i < numRot(pos); i++) {
        assig[pos] = i;
        exEn(assig, pos+1);
    }
}

// perform DEE for given sidechain in given neighbourhood. remove rotamers that are dominated by others for that sidechain
void SCenergyProvider::deeGoldstein(int sci) {
    if(numRot(sci) <= 1) return; // TODO
    // find the SCIs which inteact with sci
    vector<int> SCIs;
    for(int sci1 = 0; sci1 < chibs.size(); sci1++)
        if(sci != sci1 && interact(sci,sci1))
            SCIs.push_back(sci1);
    set<int> dominated;
    for(int ri=0; ri < numRot(sci); ri++)
        for(int rj=0; rj < numRot(sci); rj++) {
            if(ri == rj) continue;
            if(dominated.find(rj) != dominated.end()) continue;
            float Ediff = selfEn(sci, ri) - selfEn(sci, rj);
            for(vector<int>::iterator scit=SCIs.begin(); scit != SCIs.end(); ++scit) {
                float mindiff = pairEn(*scit, sci, 0, ri) - pairEn(*scit, sci, 0, rj);
                for(int rk=0; rk<numRot(*scit); rk++) {
                    float e = pairEn(*scit, sci, rk, ri) - pairEn(*scit, sci, rk, rj);
                    if(mindiff > e) mindiff = e;
                }
                Ediff += mindiff;
            }
            if(Ediff > 0) {
                dominated.insert(ri);
                //cout << ri << " -> " << rj << endl;
            }
        }
    cout << dominated.size() << " dominated of " << numRot(sci) <<" in "<< sci << endl;
    vector<int> newValid; vector<float> newEself;
    map<int,int> old2new;
    for(int i=0; i < numRot(sci); i++) {
        if(dominated.find(i) == dominated.end()) { // not dominated
            newValid.push_back( validRotamers[sci][i] );
            newEself.push_back(eself[sci][i]);
            old2new[i] = newEself.size()-1;
        }
    }
    //printEpair();
    for(int sci1=0; sci1 < sci; sci1++) {
        if(!interact(sci1,sci)) continue;
        for(int k=0; k < numRot(sci1); k++) {
            map<int,float> newepair;
            for(int i=0; i < numRot(sci); i++)
                if(old2new.find(i) != old2new.end()) // not dominated
                    newepair[ old2new[i] ] = epair[sci1][sci][k][i];
            epair[sci1][sci][k] = newepair;
        }
        calcIEpairmin(sci1,sci);
    }
    for(int sci1=sci+1; sci1 < chibs.size(); sci1++) {
        if(!interact(sci,sci1)) continue;
        map<int, map<int, float> > newepair;
        for(int i=0; i < numRot(sci); i++)
            if(old2new.find(i) != old2new.end()) // not dominated
                newepair[ old2new[i] ] = epair[sci][sci1][i];
        epair[sci][sci1] = newepair;
        calcIEpairmin(sci,sci1);
    }

    for(int i=0; i < newValid.size(); i++)
        cout << newValid[i] <<":"<< newEself[i]<<" ";
    cout << " COMPARE " << sci << endl;
    for(int i=0; i < validRotamers[sci].size(); i++)
        cout << validRotamers[sci][i] <<":"<<eself[sci][i]<<" "; cout << " COMPARE " << sci << endl;
    validRotamers[sci] = newValid;
    eself[sci] = newEself;
    calcEselfmin(sci);
}

void SCenergyProvider::printEpair() {
    cout << "------------------------------------------------------" << endl;
    for(map<int, map<int, map<int, map<int, float> > > >::iterator scit = epair.begin(); scit != epair.end(); scit++)
        for(map<int, map<int, map<int, float> > >::iterator scit1 = scit->second.begin(); scit1 != scit->second.end(); scit1++)
            for(map<int, map<int, float> >::iterator rit = scit1->second.begin(); rit != scit1->second.end(); rit++)
                for(map<int, float>::iterator rit1 = rit->second.begin(); rit1 != rit->second.end(); rit1++)
                    cout << scit->first <<" "<< scit1->first <<" "<< rit->first <<" "<< rit1->first <<" "<< rit1->second << endl;
    cout << "------------------------------------------------------" << endl;
}
//float SCenergyProvider::interCBdist(int sci1, int sci2) {
//    int CBindex1 = chibs[sci1]->getIP()[4];
//    int CBindex2 = chibs[sci2]->getIP()[4];
//    return calcDist((*pts)[CBindex1], (*pts)[CBindex2]);
//}
//float SCenergyProvider::getMaxCBDist(int sci) {
//    if(maxDistFromCB.size() == 0) numRot(0);
//    return maxDistFromCB[sci];
//}
