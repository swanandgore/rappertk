#include "EDrestraint.h"

float EDrestraint::negpen = 3.;
float EDrestraint::getPenalty() { return negpen; }
void EDrestraint::setPenalty(float np) { negpen = np; }

float EDrestraint::scatterRad = 1.;
float EDrestraint::getScatRad() { return scatterRad; }
void EDrestraint::setScatRad(float sr) { scatterRad = sr; }

//{{{ some dumb implementation for compiling without clipper
#ifndef CLIPPER
SingleMap * SingleMap::instance() { return NULL; }
void * SingleMap::map(const char *mtzfn, const char *folabel, const char *fclabel, const char *philabel, int maptype, float& mean, float& stdev) { return NULL; }
void * SingleMap::map(const char *mapfilename, float & mean, float & stdev) { return NULL; }
void SingleMap::addMap(string mtzfn, void *xmap) {}

float EDrestraint::findBadfit(const char *mtzfn, char *folabel, char *fclabel, char *philabel, vector<vector<float> > & pts, vector<int> & ptinds) {}
EDrestraint::EDrestraint(vector<int> & pis, const char *desc) : Restraint(pis,desc) {}
EDrestraint EDrestraint::makeEDrestraintFromMTZ(vector<int>& pis, const char *desc, const char *mtzfn, const char *folabel, const char *fclabel, const char *philabel,
        int maptype, float zero, float max, float cutoff)
    { return EDrestraint(pis,desc); }
EDrestraint EDrestraint::makeEDrestraintFromMap(vector<int>& pis, const char *desc, const char *mapfn, float zero, float max, float cutoff)
    { return EDrestraint(pis,desc); }
bool EDrestraint::satisfied(VVF & pts) { return false; }
float EDrestraint::score(VVF & pts) { return -999; }
void EDrestraint::describe() {}
#else
//}}}

#include "misc/verbosity.h"

#include <vector>
using std::vector;
#include <set>
using std::set;
#include "clipper/clipper-ccp4.h"
#include "clipper/clipper.h"
#include "clipper/core/clipper_types.h"
using clipper::Xmap;

struct compareCoordGrid
{
  bool operator()(const clipper::Coord_grid & s1, const clipper::Coord_grid & s2) const
  {
    if (s1.u() < s2.u()) return true;
    else if (s1.u() > s2.u()) return false;
    if (s1.v() < s2.v()) return true;
    else if (s1.v() > s2.v()) return false;
    if (s1.w() < s2.w()) return true;
    else if (s1.w() > s2.w()) return false;
    return false;
  }
};
//typedef vector<clipper::Coord_grid> SetCoordGrid;
typedef set<clipper::Coord_grid, compareCoordGrid> SetCoordGrid;

// find set of Coord_grid's within the envelope defined by given orthogonal coordinates and radius
SetCoordGrid findCoordGrids(vector<vector<float> > & pts, vector<int> & ptinds, float radius, const clipper::Cell & cell, const clipper::Grid_sampling gs) {
    const clipper::Grid_range gr(cell, gs, radius);
    SetCoordGrid scg;
    for(int pi=0; pi < ptinds.size(); pi++) {
        clipper::Coord_orth co ( pts[ptinds[pi]][0], pts[ptinds[pi]][1], pts[ptinds[pi]][2] );
        clipper::Coord_frac cf = co.coord_frac(cell);
        clipper::Coord_grid cg = cf.coord_grid(gs);
        for ( int i = gr.min().u(); i <= gr.max().u(); i++ )
            for ( int j = gr.min().v(); j <= gr.max().v(); j++ )
                for ( int k = gr.min().w(); k <= gr.max().w(); k++ ) {
                    clipper::Coord_frac tcf = (cg + clipper::Coord_grid(i,j,k)).coord_frac(gs);
                    clipper::Coord_orth tco = tcf.coord_orth(cell);
                    float dist = sqrt( (co.x()-tco.x())*(co.x()-tco.x()) + (co.y()-tco.y())*(co.y()-tco.y()) + (co.z()-tco.z())*(co.z()-tco.z()) );
                    if(dist > radius) continue;
                    scg.insert( cg + clipper::Coord_grid(i,j,k) );
                    //scg.push_back( cg + clipper::Coord_grid(i,j,k) );
                    //cout << co.x() <<" "<< co.y() <<" "<< co.z() << endl;
                    //cout << tco.x() <<" "<< tco.y() <<" "<< tco.z() <<" "<< dist <<" "<< radius << endl;
                    //if(dist > radius) exit(1);
                }
    }
    return scg;
}

bool comparableMaps(clipper::Xmap<float> & xmap1, clipper::Xmap<float> & xmap2) {
    if(xmap1.grid_sampling() != xmap2.grid_sampling() || ! xmap1.cell().equals(xmap2.cell()) ) { return false; } //TODO
    return true;
}

float corrCoeff(vector<float> & z1, vector<float> & z2) {
    float sq1=0, sq2=0, prod=0;
    for(int i=0; i < z1.size(); i++) {
        sq1 += z1[i]*z1[i];
        sq2 += z2[i]*z2[i];
        prod += z1[i]*z2[i];
    }
    return prod / sqrt(sq1*sq2);
}

float findMapCorrel(vector<vector<float> > & pts, vector<int> & ptinds, float radius,
        clipper::Xmap<float> & xmap1, float mean1, float stdev1,
        clipper::Xmap<float> & xmap2, float mean2, float stdev2) {
    if( ! comparableMaps(xmap1,xmap2) ) {
        cout << "Maps cant be compared due to different cells or grid-samplings" << endl; exit(1);
    }
    SetCoordGrid scg = findCoordGrids(pts, ptinds, radius, xmap1.cell(), xmap2.grid_sampling());
    int count = 0, cnt = 0; float rsrN = 0, rsrD = 0;
    vector<float> z1, z2;
    for(SetCoordGrid::iterator it=scg.begin(); it != scg.end(); it++) {
        float e1 = xmap1.get_data(*it);
        float e2 = xmap2.get_data(*it);
        if((e1-mean1)/stdev1 < 0.5 && (e2-mean2)/stdev2 < 0.5) continue; //both densities insignificant
        //if(e1 < 0.0001) continue; //calc-map must be nonzero
        if(fabs(e1) > 1e5 || fabs(e2) > 1e5) { cnt++; continue; } //probably NaN
        count++;
        rsrN += fabs(e1-e2); rsrD += fabs(e1+e2);


        //clipper::Coord_orth corth = it->coord_frac( xmap1.grid_sampling() ).coord_orth( xmap1.cell() );
        //clipper::Coord_orth corth1 = it->coord_frac( xmap2.grid_sampling() ).coord_orth( xmap2.cell() );
        //cout << corth[0] <<" "<< corth[1] <<" "<< corth[2] <<" : "<< corth1[0] <<" "<< corth1[1] <<" "<< corth1[2] <<" : "<< (e1-mean1)/stdev1 <<" "<< (e2-mean2)/stdev2 << endl;
        //if(e1-mean1 < 0 || e2-mean2 < 0) cout << e1 <<" "<< e2 << " "<< e1-mean1 <<" "<< e2-mean2 << endl;

        if(e1-mean1 < 0.) z1.push_back(0.); else z1.push_back( (e1-mean1)/stdev1 );
        if(e2-mean2 < 0.) z2.push_back(0.); else z2.push_back( (e2-mean2)/stdev2 );
    }
    float xcorr = corrCoeff(z1,z2); float rsr = rsrN/rsrD;
    //    if(verbose(6)) cout << "DENSITY-CORRELATION " << xcorr <<" RSR "<< rsr <<" over "<< count<<"("<<cnt<<")" <<" gridpoints for "<< ptinds.size() <<" coordinates."<< endl;
    if(count <= 0) return -1;
    return xcorr;
}

SingleMap* SingleMap::instance() {
    static SingleMap * inst = new SingleMap();
    return inst;
}

clipper::String findPrefixFor(vector<clipper::String> & allLabels, clipper::String aLabel1, clipper::String aLabel2) {
    clipper::String prefix1, prefix2;
    for(int i=0; i < allLabels.size(); i++) {
        vector<clipper::String> toks = allLabels[i].split(clipper::String(" "));
        for(int j=0; j < toks.size(); j++) {
            if(toks[j].tail() == aLabel1) prefix1 = toks[j].notail();
            if(toks[j].tail() == aLabel2) prefix2 = toks[j].notail();
        }
    }
    if(prefix1 != prefix2) { cout << "diff prefixes : " << prefix1 <<" "<< prefix2 << endl; exit(1); }
    return prefix1 + "/[" + aLabel1 + "," + aLabel2 + "]";
}

clipper::String labelFor(vector<clipper::String> & allLabels, clipper::String aLabel) {
    for(int i=0; i < allLabels.size(); i++) {
        cout << allLabels[i] << endl;
        vector<clipper::String> toks = allLabels[i].split(clipper::String("/ "));
        for(int j=0; j < toks.size(); j++) {
            cout << toks[j] << endl;
            if(toks[j] == aLabel) return allLabels[i];
        }
    }
    cout << "Label " << aLabel << " not found in mtz file. Fatal error." << endl;
    exit(1);
    return clipper::String("LABEL-NOT-FOUND");
}

void SingleMap::addMap(string mapkey, void *xmap1) {
    clipper::Xmap<float> *xmap = (clipper::Xmap<float>*) xmap1;
    xmaps[mapkey] = xmap;
    float mean = 0, stdev = 0;
    int count = 0;
    const clipper::Grid_range & gr = xmap->grid_asu();
    clipper::Coord_grid start = gr.min(), stop = gr.max();
    clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;
    i0 = clipper::Xmap_base::Map_reference_coord(*xmap,start);
    for ( iu = i0; iu.coord().u() <= stop.u(); iu.next_u() )
        for ( iv = iu; iv.coord().v() <= stop.v(); iv.next_v() )
            for ( iw = iv; iw.coord().w() <= stop.w(); iw.next_w() ) {
                mean += (*xmap)[iw]; count++;
            }
    mean /= count;
    mapMeans[mapkey] = mean;
    for ( iu = i0; iu.coord().u() <= stop.u(); iu.next_u() )
        for ( iv = iu; iv.coord().v() <= stop.v(); iv.next_v() )
            for ( iw = iv; iw.coord().w() <= stop.w(); iw.next_w() )
                stdev += ((*xmap)[iw] - mean) * ((*xmap)[iw] - mean);
    stdev = sqrt(stdev/count);
    mapStdevs[mapkey] = stdev;
}


void*  SingleMap::map(const char *mapfilename, float & mean, float & stdev) {
    string mapfn(mapfilename);
    if(xmaps.find(mapfn) == xmaps.end()) {
        clipper::Xmap<float> *xmap = new clipper::Xmap<float>();
        clipper::CCP4MAPfile file;
        file.open_read(mapfn.c_str());
        file.import_xmap( *xmap );
        file.close_read();
        addMap(mapfn, xmap);
    }
    mean = mapMeans[mapfn];
    stdev = mapStdevs[mapfn];
    return xmaps[mapfn];
}

void removeMissing(clipper::HKL_info & hklinfo, clipper::HKL_data<clipper::data32::F_phi> & hkldata) {
    vector<clipper::HKL> vechkl;
    vector<clipper::data32::F_phi> vecfphi;
    int missing = 0;
    for(clipper::HKL_info::HKL_reference_index ih = hkldata.first(); !ih.last(); ih.next()) {
        if(hkldata[ih].missing()) { missing++; continue; }
        //if(freeR[ih].flag() == 0) { rflagged++; continue; }
        vechkl.push_back(ih.hkl());
        vecfphi.push_back(hkldata[ih]);
    }
    cout << "missing " << missing <<" "<< vechkl.size() <<" "<< vecfphi.size() << endl;

    clipper::HKL_info newinfo(hklinfo.spacegroup(), hklinfo.cell(), hklinfo.resolution());
    newinfo.add_hkl_list(vechkl);
    hklinfo = newinfo;

    clipper::HKL_data<clipper::data32::F_phi> newdata(hklinfo);
    missing = 0;
    for(int i=0; i < vechkl.size(); i++) {
        if(!newdata.set_data(vechkl[i], vecfphi[i])) cout << "cant assign" << endl;
    }
    hkldata = newdata;
}

float EDrestraint::findBadfit1(const char *mapfn1, const char *mapfn2, vector<vector<float> > & pts, vector<int> & ptinds) {
    float meanFP=0, stdevFP=0, meanFC=0, stdevFC=0;
    clipper::Xmap<float> * xmapFP = (clipper::Xmap<float>*) SingleMap::instance()->map(mapfn1, meanFP, stdevFP);
    clipper::Xmap<float> * xmapFC = (clipper::Xmap<float>*) SingleMap::instance()->map(mapfn2, meanFC, stdevFC);
    return findMapCorrel(pts, ptinds, scatterRad, *xmapFP, meanFP, stdevFP, *xmapFC, meanFC, stdevFC);
}

float EDrestraint::findBadfit(const char *mtzfn, char *fplabel, char *fclabel, char *philabel, vector<vector<float> > & pts, vector<int> & ptinds) {
    float meanFP=0, stdevFP=0, meanFC=0, stdevFC=0;
    clipper::Xmap<float> * xmapFP = (clipper::Xmap<float>*) SingleMap::instance()->map(mtzfn, fplabel, fclabel, philabel, 0, meanFP, stdevFP);
    clipper::Xmap<float> * xmapFC = (clipper::Xmap<float>*) SingleMap::instance()->map(mtzfn, fclabel, fplabel, philabel, 1, meanFC, stdevFC);
    return findMapCorrel(pts, ptinds, scatterRad, *xmapFP, meanFP, stdevFP, *xmapFC, meanFC, stdevFC);
}

void*  SingleMap::map(const char *mtzfn, const char *folabel, const char *fclabel, const char *philabel, int maptype, float & mean, float & stdev) {
    char *mapstring = "ABCDEFGHIJKLMNOPQRSTUV";
    string mapkey( string(mtzfn) + folabel + fclabel + philabel + "_" + mapstring[maptype] + ".map" );
    if(xmaps.find(mapkey) == xmaps.end()) {
        clipper::CCP4MTZfile mtzin;
        clipper::HKL_info myhkl;
        clipper::HKL_data<clipper::data32::F_phi> fphidata(myhkl), fcdata(myhkl);
        mtzin.open_read(mtzfn);
        mtzin.import_hkl_info(myhkl); // read spacegroup, cell, resolution, hkls
        vector<clipper::String> clabels = mtzin.column_labels();
        clipper::String foPhiLabel = findPrefixFor(clabels, clipper::String(folabel), clipper::String(philabel));
        clipper::String fcPhiLabel = findPrefixFor(clabels, clipper::String(fclabel), clipper::String(philabel));
        cout << "F1-PHIC = " << foPhiLabel << endl;
        mtzin.import_hkl_data(fphidata, foPhiLabel);
        cout << "F2-PHIC = " << fcPhiLabel << endl;
        mtzin.import_hkl_data(fcdata  , fcPhiLabel);
        mtzin.close_read();
        if(maptype == 0) { // 2 Fo - Fc
            fphidata = fphidata + fphidata - fcdata;
            cout << "2F1-F2 map" << endl;
        }
        else if(maptype == 1) { // Fo
            cout << "F1 map" << endl;
        }
        else { cout << "Unknown maptype " << maptype << endl; exit(0); }
        cout << "NUMOBS " << myhkl.num_reflections() <<" "<< fphidata.num_obs() << endl;
        removeMissing(myhkl, fphidata);
        cout << "NUMOBS " << myhkl.num_reflections() <<" "<< fphidata.num_obs() << endl;
        clipper::Resolution resol( 1.0/sqrt(fphidata.invresolsq_range().max()) );
        clipper::Grid_sampling mygrid( myhkl.spacegroup(), myhkl.cell(), resol );  // define grid
        clipper::Xmap<float> *xmap = new clipper::Xmap<float> ( myhkl.spacegroup(), myhkl.cell(), mygrid );  // define map
        cout << "Grid..." << xmap->grid_sampling().format() << "\n";
        xmap->fft_from(fphidata); // generate map
        addMap(mapkey, xmap);

        clipper::CCP4MAPfile myfile;
        myfile.open_write( mapkey.c_str() ); //"rtk.map" );
        myfile.export_xmap( *xmap );
        myfile.close_write();
        cout << "Wrote map to " << mapkey << endl;

    }
    mean = mapMeans[mapkey];
    stdev = mapStdevs[mapkey];
    return xmaps[mapkey];
}

EDrestraint::EDrestraint(vector<int> & pis, const char *desc) : Restraint(pis,desc) {}

EDrestraint EDrestraint::makeEDrestraintFromMTZ(vector<int>& pis, const char *desc, const char *mtzfn, const char *folabel, const char *fclabel, const char *philabel,
        int maptype, float zeroDenSig, float maxDenSig, float meanDenSig) {
    EDrestraint edr(pis, desc);
    edr.xmap = SingleMap::instance()->map(mtzfn, folabel, fclabel, philabel, maptype, edr.mean, edr.stdev); 
    edr.zero = edr.mean + zeroDenSig * edr.stdev;
    edr.max = edr.mean + maxDenSig * edr.stdev;
    edr.cutoff = edr.mean + meanDenSig * edr.stdev;
    //    if(verbose(7)) cout << mtzfn<<" "<<zeroDenSig<<" "<<maxDenSig<<" "<<meanDenSig<<" mean,stdev "<< edr.mean <<" "<< edr.stdev <<" "<< edr.zero <<" "<< edr.max <<" "<< edr.cutoff<< endl;
    return edr;
}

EDrestraint EDrestraint::makeEDrestraintFromMap(vector<int>& pis, const char *desc, const char *mapfn, float zeroDenSig, float maxDenSig, float meanDenSig) {
    EDrestraint edr(pis, desc);
    edr.xmap = SingleMap::instance()->map(mapfn, edr.mean, edr.stdev); 
    edr.zero = edr.mean + zeroDenSig * edr.stdev;
    edr.max = edr.mean + maxDenSig * edr.stdev;
    edr.cutoff = edr.mean + meanDenSig * edr.stdev;
    //    if(verbose(7)) cout << mapfn<<" "<<zeroDenSig<<" "<<maxDenSig<<" "<<meanDenSig<<" mean,stdev "<< edr.mean <<" "<< edr.stdev <<" "<< edr.zero <<" "<< edr.max <<" "<< edr.cutoff<< endl;
    return edr;
}

float EDrestraint::scoreSigma(VVF & pts) {
    return (score(pts) - mean) / stdev;
}

float EDrestraint::evalDensity(float den) {
    if(den > max) return max;
    if(den < zero) return (-1 * negpen * (zero - den) + zero);
    return den;
}

//#define TESTDEVELOP 1
#ifdef TESTDEVELOP
float EDrestraint::scoreAround(VVF & pts) {
    float scatrad = 2.;
    SetCoordGrid scg = findCoordGrids( pts, ptInds, scatrad,
        ((clipper::Xmap<float>*)xmap)->cell(),
        ((clipper::Xmap<float>*)xmap)->grid_sampling() );
    vector<float> mapdensity, shape;
    int count = -1, pi=0;
    for(SetCoordGrid::iterator it=scg.begin(); it != scg.end(); it++) {
        float mapdata = ((clipper::Xmap<float>*)xmap)->get_data(*it);
        if( ((clipper::Xmap<float>*)xmap)->get_data(*it) < 0.5 ) continue;
        if( ((clipper::Xmap<float>*)xmap)->get_data(*it) > 10 ) continue;
        mapdensity.push_back( ((clipper::Xmap<float>*)xmap)->get_data(*it) );
        shape.push_back(0.);
        count ++;
        clipper::Coord_orth corth = it->coord_frac( ((clipper::Xmap<float>*)xmap)->grid_sampling() ).coord_orth( ((clipper::Xmap<float>*)xmap)->cell() );
        for(int i=0; i < ptInds.size(); i++) {
            pi = ptInds[i];
            float dist = sqrt((pts[pi][0]-corth[0])*(pts[pi][0]-corth[0]) + (pts[pi][1]-corth[1])*(pts[pi][1]-corth[1]) + (pts[pi][2]-corth[2])*(pts[pi][2]-corth[2]));
            if(dist <= scatrad) shape[count] += scatrad-dist;
        }
    }
    //cout << corrCoeff(mapdensity, shape) <<" "<< mapdensity.size() <<" "<< ptInds.size() << endl;
    if(ptInds.size() * 4 > mapdensity.size()) return 0;
    return corrCoeff(mapdensity, shape);
}
#else//TESTDEVELOP
float EDrestraint::scoreAround(VVF & pts) {
    SetCoordGrid scg = findCoordGrids( pts, ptInds, scatterRad,
        ((clipper::Xmap<float>*)xmap)->cell(),
        ((clipper::Xmap<float>*)xmap)->grid_sampling() );
    float tot = 0.;
    for(SetCoordGrid::iterator it=scg.begin(); it != scg.end(); it++) {
        //clipper::Coord_orth corth = it->coord_frac( ((clipper::Xmap<float>*)xmap)->grid_sampling() ).coord_orth( ((clipper::Xmap<float>*)xmap)->cell() );
        //cout << corth[0] <<" "<< corth[1] <<" "<< corth[2] <<" ";
        //cout << ((clipper::Xmap<float>*)xmap)->get_data(*it) <<" "<< evalDensity( ((clipper::Xmap<float>*)xmap)->get_data(*it) ) << endl;
        tot += evalDensity( ((clipper::Xmap<float>*)xmap)->get_data(*it) );
    }
    //cout << "EVALDEN " << tot << endl;
    return tot;
}
#endif//TESTDEVELOP

float EDrestraint::score(VVF & pts) {
    float avgSigma = 0, density = 0;
    for(int i=0; i < ptInds.size(); i++) {
        clipper::Coord_orth corth(pts[ptInds[i]][0], pts[ptInds[i]][1], pts[ptInds[i]][2]);
        clipper::Coord_frac cfrac = corth.coord_frac(((clipper::Xmap<float>*)xmap)->cell());
        clipper::Coord_grid cgrid = cfrac.coord_grid(((clipper::Xmap<float>*)xmap)->grid_sampling());
        density = ((clipper::Xmap<float>*)xmap)->get_data(cgrid);
	//        if(verbose(7)) cout << ptInds[i] << " density " << density;
        density = evalDensity(density);
	//        if(verbose(7)) cout << " " << density << endl;
        avgSigma += density;
    }
    //    if(verbose(7)) cout << "SUMSIG " << avgSigma <<" "<< ptInds.size() <<" "<< avgSigma/ptInds.size() << endl;
    return avgSigma;
}

bool EDrestraint::satisfied(VVF & pts) {
  //    if(verbose(7)) cout << "checking EDrestraint " << zero <<" "<< max <<" "<< mean <<" "<< stdev << endl;
    float avgSigma = score(pts) / ptInds.size();
    if(avgSigma < cutoff) return false;
    //    if(verbose(7)) cout << "true" << endl;
    return true;
}

void EDrestraint::describe() {
    cout << "EDrestraint" << endl; 
}

#endif//CLIPPER
