// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <unordered_map>

#include "Interface.h"

using namespace std;

Interface::Interface(Run * run) {
    run_ = run;
    s0_ = 5.40; // 0~5.82
    energy_ = 0.;
}

int     Interface::updateForces() {
    // reset all interfaceForce values in vertices
    for (long int i = 0; i < run_->vertices_.size(); i++) {
        for (int j = 0; j < 3; j++) {
            run_->vertices_[i]->interfaceForce_[j] = 0.;
        }
    }

    // update area of each polygon
    for (auto polygon : run_->polygons_) {
        polygon->updateArea();
    }

    // update tension in each polygon
    updateTension();

    // update interfaceForce values
    for (long int i = 0; i < run_->polygons_.size(); i++) {
        updatePolygonForces(run_->polygons_[i]);
    }

    return 0;
}

int Interface::updatePolygonForces(Polygon *polygon) {
    // reset interfaceForce values of the polygon center
    for (int m = 0; m < 3; m++) {
        polygon->interfaceForce_[m] = 0.;
    }
    double tension = polygon->tension_;

    // the polygon center is the reference point
    for (int i = 0; i < polygon->edges_.size(); i++) {
        Edge * edge = polygon->edges_[i];
        // the vectors pointing from polygon center to edge vertices
        double cv[2][3];
        for (int k = 0; k < 2; k++) {
            Vertex * vertex = edge->vertices_[k];
            for (int m = 0; m < 3; m++) {
                cv[k][m] = vertex->position_[m] - polygon->center_[m];
            }
            while (cv[k][0] > run_->Lx_/2.0) {
                cv[k][0] = cv[k][0] - run_->Lx_;
            }
            while (cv[k][0] < (-1.0)*run_->Lx_/2.0) {
                cv[k][0] = cv[k][0] + run_->Lx_;
            }
            while (cv[k][1] > run_->Ly_/2.0) {
                cv[k][1] = cv[k][1] - run_->Ly_;
            }
            while (cv[k][1] < (-1.0)*run_->Ly_/2.0) {
                cv[k][1] = cv[k][1] + run_->Ly_;
            }
            while (cv[k][2] > run_->Lz_/2.0) {
                cv[k][2] = cv[k][2] - run_->Lz_;
            }
            while (cv[k][2] < (-1.0)*run_->Lz_/2.0) {
                cv[k][2] = cv[k][2] + run_->Lz_;
            }
        }
        // the edge vector
        double vv[3];
        for (int m = 0; m < 3; m++) {
            vv[m] = cv[1][m] - cv[0][m];
            // vv[m] = edge->vv_[m];
        }
        // compute the normal vector of the triangle interface formed by polygon center, and edge vertices
        double nv[3];
        nv[0] = cv[0][1]*cv[1][2] - cv[1][1]*cv[0][2];
        nv[1] = cv[1][0]*cv[0][2] - cv[0][0]*cv[1][2];
        nv[2] = cv[0][0]*cv[1][1] - cv[1][0]*cv[0][1];
        double norm_nv = sqrt(nv[0]*nv[0] + nv[1]*nv[1] + nv[2]*nv[2]);
        nv[0] = nv[0] / norm_nv;
        nv[1] = nv[1] / norm_nv;
        nv[2] = nv[2] / norm_nv;
        // compute forces on triangle edges
        double Fcv0[3];
        double Fcv1[3];
        double Fvv[3];
        Fvv[0] = tension*(nv[1]*vv[2] - vv[1]*nv[2]);
        Fvv[1] = tension*(vv[0]*nv[2] - nv[0]*vv[2]);
        Fvv[2] = tension*(nv[0]*vv[1] - vv[0]*nv[1]);
        Fcv0[0] = tension*(nv[1]*cv[0][2] - cv[0][1]*nv[2]);
        Fcv0[1] = tension*(cv[0][0]*nv[2] - nv[0]*cv[0][2]);
        Fcv0[2] = tension*(nv[0]*cv[0][1] - cv[0][0]*nv[1]);
        Fcv1[0] = tension*(cv[1][1]*nv[2] - nv[1]*cv[1][2]);
        Fcv1[1] = tension*(nv[0]*cv[1][2] - cv[1][0]*nv[2]);
        Fcv1[2] = tension*(cv[1][0]*nv[1] - nv[0]*cv[1][1]);
        // update interfaceForces
        for (int m = 0; m < 3; m++) {
            edge->vertices_[0]->interfaceForce_[m] = edge->vertices_[0]->interfaceForce_[m] + 0.5*(Fcv0[m]+Fvv[m]);
            edge->vertices_[1]->interfaceForce_[m] = edge->vertices_[1]->interfaceForce_[m] + 0.5*(Fcv1[m]+Fvv[m]);
            polygon->interfaceForce_[m] = polygon->interfaceForce_[m] + 0.5*(Fcv0[m]+Fcv1[m]);
        }
    }

    // redistribute polygon center interfaceForces back to vertices
    double sum_l = 0.;
    for (int i = 0; i < polygon->edges_.size(); i++) {
        sum_l += polygon->edges_[i]->length_;
    }
    for (int i = 0; i < polygon->edges_.size(); i++) {
        double weight = polygon->edges_[i]->length_/sum_l;
        for (int k = 0; k < 2; k++) {
            Vertex *vertex = polygon->edges_[i]->vertices_[k];
            for (int m = 0; m < 3; m++) {
                vertex->interfaceForce_[m] = vertex->interfaceForce_[m] + 0.5*weight*polygon->interfaceForce_[m];
            }
        }
    }

    return 0;
}

int Interface::updateTension() {
    for (auto polygon : run_->polygons_) {
        polygon->tension_ = 0.;
    }
    for (auto cell : run_->cells_) {
        double s = 0.;
        int springtrigger = 0;
        for (auto polygon : cell->polygons_) {
            s += polygon->area_;
            if(cell->type_==31){springtrigger = 1;}
        }
        for (auto polygon : cell->polygons_) {
            double effsigma = 0.;
            //Basement = 1 (0), Basak = 2 (1), Suprabasal = 3 (3)
            //FinalBasal = 4 (7), FinalSupra = 5 (15), FinalBase = 6 (31), Delam Cell = 7 (63)
            if(polygon->type_==1){effsigma = run_->sigma12_;}
            if(polygon->type_==3){effsigma = run_->sigma31_;}
            if(polygon->type_==4){effsigma = run_->sigma23_;}
            if(polygon->type_==7){effsigma = run_->sigma14_;}
            if(polygon->type_==8){effsigma = run_->sigma24_;}
            if(polygon->type_==10){effsigma = run_->sigma34_;}
            if(polygon->type_==15){effsigma = run_->sigma15_;}
            if(polygon->type_==16){effsigma = run_->sigma25_;}
            if(polygon->type_==18){effsigma = run_->sigma35_;}
            if(polygon->type_==22){effsigma = run_->sigma45_;}
            if(polygon->type_==31){effsigma = run_->sigma16_;}
            if(polygon->type_==32){effsigma = run_->sigma26_;}
            if(polygon->type_==34){effsigma = run_->sigma36_;}
            if(polygon->type_==38){effsigma = run_->sigma46_;}
            if(polygon->type_==46){effsigma = run_->sigma56_;}
            if(polygon->type_==63){effsigma = run_->sigma17_;}
            if(polygon->type_==64){effsigma = run_->sigma27_;}
            if(polygon->type_==66){effsigma = run_->sigma37_;}
            if(polygon->type_==70){effsigma = run_->sigma47_;}
            if(polygon->type_==78){effsigma = run_->sigma57_;}
            if(polygon->type_==94){effsigma = run_->sigma67_;}

            //Tension from interface
            if(springtrigger==0){polygon->tension_ += 2.0*(s - cell->s0_)+effsigma;}
            else{polygon->tension_ += 2.0*(s - cell->s0_)+effsigma;}
            
        }
    }


    return 0;
}

int Interface::updateEnergy() {
    energy_ = 0.;
    int counts = 0.;

    for (auto cell : run_->cells_) {
        for (auto polygon : cell->polygons_) {
            polygon->type_ = 0;
        }
    }

    for (auto cell : run_->cells_) {
        for (auto polygon : cell->polygons_) {
            polygon->type_ += cell->type_;
        }
    }

    for (auto cell : run_->cells_) {
        double s = 0.;
        double HITE = 0.;
        int springtrigger = 0;

        if(cell->type_==31){springtrigger = 1;}

        for (auto polygon : cell->polygons_) {
            s += polygon->area_;
            //Basement = 1 (0), Basal = 2 (1), Suprabasal = 3 (3)
            //FinalBasal = 4 (7), FinalSupra = 5 (15), FinalBase = 6 (31), Delam Cell = 7 (63)
            if(polygon->type_==1){HITE += run_->sigma12_*polygon->area_;}
            if(polygon->type_==3){HITE += run_->sigma31_*polygon->area_;} 
            if(polygon->type_==4){HITE += run_->sigma23_*polygon->area_;}
            if(polygon->type_==7){HITE += run_->sigma14_*polygon->area_;}
            if(polygon->type_==8){HITE += run_->sigma24_*polygon->area_;} 
            if(polygon->type_==10){HITE += run_->sigma34_*polygon->area_;}
            if(polygon->type_==15){HITE += run_->sigma15_*polygon->area_;}
            if(polygon->type_==16){HITE += run_->sigma25_*polygon->area_;}
            if(polygon->type_==18){HITE += run_->sigma35_*polygon->area_;}
            if(polygon->type_==22){HITE += run_->sigma45_*polygon->area_;}
            if(polygon->type_==31){HITE += run_->sigma16_*polygon->area_;}
            if(polygon->type_==32){HITE += run_->sigma26_*polygon->area_;}
            if(polygon->type_==34){HITE += run_->sigma36_*polygon->area_;}
            if(polygon->type_==38){HITE += run_->sigma46_*polygon->area_;}
            if(polygon->type_==46){HITE += run_->sigma56_*polygon->area_;}
            if(polygon->type_==63){HITE += run_->sigma17_*polygon->area_;}
            if(polygon->type_==64){HITE += run_->sigma27_*polygon->area_;}
            if(polygon->type_==66){HITE += run_->sigma37_*polygon->area_;}
            if(polygon->type_==70){HITE += run_->sigma47_*polygon->area_;}
            if(polygon->type_==78){HITE += run_->sigma57_*polygon->area_;}
            if(polygon->type_==94){HITE += run_->sigma67_*polygon->area_;}
        }

        //Energy with or w/o springs
        if(springtrigger==0){
            energy_ += pow(s - cell->s0_, 2.0) + HITE/2;
        }
        else{
            energy_ += pow(s - cell->s0_, 2.0) + HITE/2;
        }
        
    }

    return 0;
}