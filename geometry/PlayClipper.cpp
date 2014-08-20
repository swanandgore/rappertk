#ifdef CLIPPERGEOM

#include "PlayClipper.h"

#include <vector>
#include <iostream>
using namespace std;

#include "clipper/core/coords.h"
using clipper::Coord_orth;
using clipper::Coord_frac;

#include "clipper/core/cell.h"

#include "clipper/core/spacegroup.h"

void generateSymmpos(clipper::Cell * cell, clipper::Spacegroup * spgr,
        vector<float> & point, vector<vector<float> > & images,
        vector<float> & llb, vector<float> & rtf, float specialCutoff, float farCutoff);


PlayClipper::PlayClipper(const char *sg, float a, float b, float c, float A, float B, float C, float x, float y, float z) {
    clipper::Spacegroup * spgr = new clipper::Spacegroup( clipper::Spgr_descr(sg) );
    clipper::Cell * cell = new clipper::Cell( clipper::Cell_descr( a, b, c, A, B, C ) );
    vector<float> target;
    target.push_back(x);
    target.push_back(y);
    target.push_back(z);

    vector<vector<float> > images;
    vector<float> llb(3), rtf(3);
    for(int i=0; i < 3; i++) { llb[i] = -1e10; rtf[i] = 1e10; }

    //generateSymmpos(cell, spgr, target, images, 0.2, 50);
    generateSymmpos(cell, spgr, target, images, llb, rtf, 0.2, 50);

    images.push_back(target);
    for(int i=0; i < images.size(); i++)
        //cout << "cgo += [SPHERE, "<< images[i][0] <<","<< images[i][1] <<","<< images[i][2] <<",1,]"<< endl;
        printf("ATOM    302  CA  CYS A  66    %8.3f%8.3f%8.3f\n", images[i][0], images[i][1], images[i][2]);
}

//    Coord_frac tcf = target.coord_frac( cell );
//    Coord_frac floors;
//    for(int i=0; i < 3; i++) floors[i] = floor(tcf[i]);
//    for(int i=0; i < 3; i++) tcf[i] -= floors[i];
//    for ( int isym = 0; isym < spgr.num_symops(); isym++ ) {
//        Coord_frac scf = spgr.symop(isym) * tcf;
//        for(int i=0; i < 3; i++) scf[i] += floors[i];
//        Coord_orth image = scf.coord_orth(cell);
//        cout << "("<< scf[0] <<" "<< scf[1] <<" "<< scf[2] <<") ("<< image[0] <<" "<< image[1] <<" "<< image[2] <<")"<< endl;
//        VecFloat im1(3); im.push_back(image[0]); im.push_back(image[1]); im.push_back(image[2]);
//        for(int i=-2; i < 3; i++) for(int j=-2; j < 3; j++) for(int k=-2; k < 3; k++) {
//            im1[0] = image[0] + i*cell.a();
//            im1[1] = image[1] + j*cell.b();
//            im1[2] = image[2] + k*cell.c();
//            if( fabs(image[0]-target[0]) < 0.2 && fabs(image[1]-target[1]) < 0.2 && fabs(image[2]-target[2]) < 0.2 ) continue; // special position
//            if( fabs(im1[0] - target[0]) > 50 || fabs(im1[1] - target[1]) > 50 || fabs(im1[2] - target[2]) > 50 ) continue; // too far away
//            results.push_back(im1);
//        }
//        cout << "("<< scf[0] <<" "<< scf[1] <<" "<< scf[2] <<") ("<< image[0] <<" "<< image[1] <<" "<< image[2] <<")"<< endl;
//    }


//std::vector<Coord_orth> results;
//
//Coord_frac tcf = target.coord_frac( cell );
//for ( int isym = 0; isym < spgr.num_symops(); isym++ ) {
//   Coord_frac scf = spgr.symop(isym) * tcf;
//   scf = scf.lattice_copy_near( tcf );
//   for ( double w = -2.0; w < 2.5; w+=1.0 )
//     for ( double v = -2.0; v < 2.5; v+=1.0 )
//       for ( double u = -2.0; u < 2.5; u+=1.0 ) {
//         Coord_frac trn = Coord_frac(u,v,w) + scf;
//         if ( (trn-target).lengthsq( cell ) < dist*dist )
//           results.push_back( trn.coord_orth( cell ) );
//       }
//}
/*
//void generateSymmpos(Cell & cell, Spacegroup & spgr, vector<float> & point, vector<vector<float> > & images, float specialCutoff, float farCutoff) {
//    Coord_orth target(point[0], point[1], point[2]);
//    Coord_frac tcf = target.coord_frac( cell );
//
//    Coord_frac floors;
//    for(int i=0; i < 3; i++) floors[i] = floor(tcf[i]);
//    for(int i=0; i < 3; i++) tcf[i] -= floors[i];
//
//    for ( int isym = 1; isym < spgr.num_symops(); isym++ ) {
//        Coord_frac scf = spgr.symop(isym) * tcf;
//        for(int i=0; i < 3; i++) scf[i] += floors[i];
//        Coord_orth image = scf.coord_orth(cell);
//        cout << "first ("<< scf[0] <<" "<< scf[1] <<" "<< scf[2] <<") ("<< image[0] <<" "<< image[1] <<" "<< image[2] <<")"<< endl;
//        vector<float> im1(3);
//        for(int i=-2; i < 3; i++) for(int j=-2; j < 3; j++) for(int k=-2; k < 3; k++) {
//            im1[0] = image[0] + i*cell.a();
//            im1[1] = image[1] + j*cell.b();
//            im1[2] = image[2] + k*cell.c();
//            if( fabs(image[0]-target[0]) < specialCutoff && fabs(image[1]-target[1]) < specialCutoff && fabs(image[2]-target[2]) < specialCutoff ) continue; // special position
//            if( fabs(im1[0] - target[0]) > farCutoff || fabs(im1[1] - target[1]) > farCutoff || fabs(im1[2] - target[2]) > farCutoff ) continue; // too far away
//            images.push_back(im1);
//        }
//    }
//}
*/

#endif//CLIPPERGEOM
