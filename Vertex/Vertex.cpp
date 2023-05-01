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

#include "Vertex.h"

using namespace std;

Vertex::Vertex(Run * run, long int id) {
    run_ = run;
    id_ = id;
    for (int i = 0; i < 3; i++) {
        position_[i] = 0.;
        volumeForce_[i] = 0.;
        interfaceForce_[i] = 0.;
        velocity_[i] = 0.;
        motility_[i] = 0.;
        springforces_[i] = 0.;
    }
}

int Vertex::logCells(std::string name) {
    printf("%s %ld\n",name.c_str(), id_);
    printf("%ld",cells_.size());
    for (auto cell : cells_) {
        printf(" %ld",cell->id_);
    }
    printf("\n");

    return 0;
}

int Vertex::logEdges(std::string name) {
    printf("%s %ld\n",name.c_str(), id_);
    printf("%ld",edges_.size());
    for (auto edge : edges_) {
        printf(" %ld",edge->id_);
    }
    printf("\n");

    return 0;
}

int Vertex::updateSP(double temperature, double temperaturebot, double currtime, int placodeon) {
    double directorx = 0;
    double directory = 0;
    double directorz = 0;
    double curtemp = 0;
    int bmtrigger = 0;
    int finbastrigger = 0;
    double tempplacx = 0;
    double tempplacy = 0;
    double tempspeed = 0;
    double tempmag = 0 ;


    for (auto cell : cells_) {
        if(cell->type_==31){curtemp = temperature;}
        else{curtemp = temperature;}


        //Turn on basement springs
        if(cell->type_==31){bmtrigger=1;}
        if(cell->type_==7||cell->type_==63){finbastrigger=1;}

        directorx += curtemp*cos(cell->cellDirectors_[1])*sin(cell->cellDirectors_[0]);
        directory += curtemp*sin(cell->cellDirectors_[1])*sin(cell->cellDirectors_[0]);
        directorz += curtemp*cos(cell->cellDirectors_[0]);
        
    }

    motility_[0] = directorx/cells_.size();
    motility_[1] = directory/cells_.size();
    motility_[2] = directorz/cells_.size();

    // if(currtime > 200){
    //     //double timeadjust = (currtime-200)/40;
    //     double vraise=0.005;
    //     for (auto cell : cells_) {
    //         if(cell->type_==63){
    //         //motility_[2] = temperaturebot*timeadjust;
    //         motility_[2] = -velocity_[2]+vraise;
    //         }
    //     }
    // }

    int pushcelltrigger = 1;
    if(currtime > 75 && pushcelltrigger==1){
    	int delam = 0;
   		int basal = 0;
    	int bm = 0;
        int suprabasal = 0;

    	//double vraise=0.01;
    	//double vraise=1e-4;
        double vraise=run_->temperaturebot_;
        for (auto cell : cells_) {

        	// for (auto polygon : cell->polygons_) {
        	// 	if(polygon->type_==70||polygon->type_==78||polygon->type_==94){
        	// 		velocity_[2] = 0;
        	// 		motility_[2] = vraise;
        	// 	}
        	// }

            if(cell->type_==63){
            delam = 1;  
            //velocity_[2] = vraise; 
            if(int(currtime/0.005)%100!=0)
                {
                velocity_[2] = 0;
                motility_[2] = vraise/0.01;
                //motility_[2] = vraise;
                }

            }
            if(cell->type_==31){
                bm = 1;
            }
            if(cell->type_==7){
                basal = 1;
            }
            if(cell->type_==15){
                suprabasal = 1;
            }
        }
        if(bm==1 && suprabasal==1){
            for (int m = 0; m < 3; m++) {
                motility_[m] = 0;
                velocity_[m] = 0;
            }
        }
        // if(delam==1){
        // 	//motility_[2] = vraise;
        // 	velocity_[2] = vraise;
        // }
    }


    // if(placodeon == 1 && currtime > 200){
    //     for (auto cell : cells_) {
    //         if(cell->type_==7){

    //         tempplacx = run_->placodex_ - cell->cellcenter_[0];
    //         tempplacy = run_->placodey_ - cell->cellcenter_[1];

    //         if(tempplacx > 3*run_->Lx_/4){
    //             tempplacx = tempplacx - run_->Lx_;
    //         }
    //         if(tempplacx < -3*run_->Lx_/4){
    //             tempplacx = tempplacx + run_->Lx_;
    //         }            
    //         if(tempplacy > 3*run_->Ly_/4){
    //             tempplacy = tempplacx - run_->Ly_;
    //         }
    //         if(tempplacy < -3*run_->Ly_/4){
    //             tempplacy = tempplacy + run_->Ly_;
    //         }      

    //         tempmag = pow(pow(tempplacx,2)+pow(tempplacy,2),0.5);
    //         tempspeed = temperaturebot*exp(-(tempmag-pow(run_->placodesize_/M_PI,0.5)*pow(cell->v0_,1/3))/(2*pow(run_->placodesize_/M_PI,0.5)*0.5*pow(cell->v0_,1/3)));

    //         motility_[0] = tempspeed*tempplacx/tempmag;
    //         motility_[1] = tempspeed*tempplacy/tempmag;
    //         motility_[2] = 0;
                           
    //         }
    //     }
    // }



    return 0;
}

int Vertex::updateVpos(double currtime){

    double vraise=1e-6;
    double oldz = 0;
    int delam = 0;
    int basal = 0;
    int bm = 0;

    for (int m = 0; m < 3; m++) {
        if(m==2){
            oldz = position_[m];
        }
        position_[m] = position_[m] + velocity_[m] * run_->dt_ + run_->dt_ * motility_[m];
    }
      

    if(currtime > 125){
        for (auto cell : cells_) {
            if(cell->type_==63){
            //position_[2] = oldz+vraise*run_->dt_;
            position_[2] = oldz+vraise;
            delam = 1;   
            }
            if(cell->type_==31){
                bm = 1;
            }
            if(cell->type_==7){
                basal = 1;
            }
        }
        if(bm==1 && basal==0 && delam==0){
            for (int m = 0; m < 3; m++) {
                position_[m] = position_[m] - velocity_[m] * run_->dt_ - run_->dt_ * motility_[m];
            }
        }

    }

    return 0;
}

int Vertex::resetTensions(){

    for (int m = 0; m < 3; m++) {
        springforces_[m] = 0.;
    }


    return 0;
}

int Vertex::magForce(){

    double totalforce = 0.0;
    for (int m = 0; m < 3; m++) {
        totalforce += (velocity_[m])*(velocity_[m]);
    }
    netforce_ = pow(totalforce,0.5);
}