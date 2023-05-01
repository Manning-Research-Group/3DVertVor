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

#include "Polygon.h"

using namespace std;

Polygon::Polygon(Run * run, long int id) {
    run_ = run;
    id_ = id;
    for (int i = 0; i < 3; i++) {
        center_[i] = 0.;
        volumeForce_[i]  = 0.;
        springtension_[i] = 0.;
        interfaceForce_[i]  = 0.;
    }
    area_ = 0.;
    tension_ = 0.;
}

int Polygon::updateVertices() {
    vertices_.clear();
    std::vector<Edge *> tmp_edges = edges_;
    Vertex * currentVertex = NULL;
    for (int i = edges_.size()-1; i > 0; i--) {
        if (i == edges_.size()-1) {
            vertices_.push_back(tmp_edges[tmp_edges.size()-1]->vertices_[0]);
            currentVertex = tmp_edges[tmp_edges.size()-1]->vertices_[1];
            vertices_.push_back(currentVertex);
            tmp_edges.pop_back();
            continue;
        }
        for (int j = 0; j < tmp_edges.size(); j++) {
            if (currentVertex == tmp_edges[j]->vertices_[0]) {
                currentVertex = tmp_edges[j]->vertices_[1];
                vertices_.push_back(currentVertex);
                if (j < tmp_edges.size() - 1) {
                    tmp_edges[j] = tmp_edges[tmp_edges.size()-1];
                }
                tmp_edges.pop_back();
                break;
            }
            if (currentVertex == tmp_edges[j]->vertices_[1]) {
                currentVertex = tmp_edges[j]->vertices_[0];
                vertices_.push_back(currentVertex);
                if (j < tmp_edges.size() - 1) {
                    tmp_edges[j] = tmp_edges[tmp_edges.size()-1];
                }
                tmp_edges.pop_back();
                break;
            }
        }
    }

    return 0;
}

int Polygon::updateCenter() {
    // set reference point
    double tmp_origin[3];
    for (int i = 0; i < 3; i++) {
        tmp_origin[i] = edges_[0]->vertices_[0]->position_[i];
    }

    double sum_lx = 0.;
    double sum_ly = 0.;
    double sum_lz = 0.;
    double sum_l = 0.;
    for (int i = 0; i < edges_.size(); i++) {
        double length = edges_[i]->length_;
        double dx = edges_[i]->center_[0] - tmp_origin[0];
        double dy = edges_[i]->center_[1] - tmp_origin[1];
        double dz = edges_[i]->center_[2] - tmp_origin[2];
        while (dx > run_->Lx_/2.0) {
            dx -= run_->Lx_;
        }
        while (dx < (-1.0)*run_->Lx_/2.0) {
            dx += run_->Lx_;
        }
        while (dy > run_->Ly_/2.0) {
            dy -= run_->Ly_;
        }
        while (dy < (-1.0)*run_->Ly_/2.0) {
            dy += run_->Ly_;
        }
        while (dz > run_->Lz_/2.0) {
            dz -= run_->Lz_;
        }
        while (dz < (-1.0)*run_->Lz_/2.0) {
            dz += run_->Lz_;
        }
        sum_lx += length*dx;
        sum_ly += length*dy;
        sum_lz += length*dz;
        sum_l += length;
    }
    center_[0] = sum_lx/sum_l + tmp_origin[0];
    center_[1] = sum_ly/sum_l + tmp_origin[1];
    center_[2] = sum_lz/sum_l + tmp_origin[2];

    return 0;
}

int Polygon::updateArea() {
    // reset polygon area
    area_ = 0.;

    // the polygon center is the reference point
    for (auto edge : edges_) {
        // the vectors pointing from polygon center to edge vertices
        double cv[2][3];
        for (int k = 0; k < 2; k++) {
            Vertex *vertex = edge->vertices_[k];
            for (int m = 0; m < 3; m++) {
                cv[k][m] = vertex->position_[m] - center_[m];
            }
            while (cv[k][0] > run_->Lx_ / 2.0) {
                cv[k][0] = cv[k][0] - run_->Lx_;
            }
            while (cv[k][0] < (-1.0) * run_->Lx_ / 2.0) {
                cv[k][0] = cv[k][0] + run_->Lx_;
            }
            while (cv[k][1] > run_->Ly_ / 2.0) {
                cv[k][1] = cv[k][1] - run_->Ly_;
            }
            while (cv[k][1] < (-1.0) * run_->Ly_ / 2.0) {
                cv[k][1] = cv[k][1] + run_->Ly_;
            }
            while (cv[k][2] > run_->Lz_ / 2.0) {
                cv[k][2] = cv[k][2] - run_->Lz_;
            }
            while (cv[k][2] < (-1.0) * run_->Lz_ / 2.0) {
                cv[k][2] = cv[k][2] + run_->Lz_;
            }
        }
        // compute the normal vector of the triangle interface formed by polygon center, and edge vertices
        double nv[3];
        nv[0] = cv[0][1] * cv[1][2] - cv[1][1] * cv[0][2];
        nv[1] = cv[1][0] * cv[0][2] - cv[0][0] * cv[1][2];
        nv[2] = cv[0][0] * cv[1][1] - cv[1][0] * cv[0][1];
        double norm_nv = sqrt(nv[0] * nv[0] + nv[1] * nv[1] + nv[2] * nv[2]);
        area_ += 0.5 * norm_nv;
    }

    return 0;
}

bool Polygon::crossBoundary() {
    for (int i = 0; i < edges_.size(); i++) {
        if (edges_[i]->crossBoundary()) {
            return true;
        }
    }

    return false;
}

bool Polygon::checkH() {
    if (edges_.size() != 3) {
        return false;
    }
    for (int j = 0; j < 3; j++) {
        if (!edges_[j]->checkH()) {
            return false;
        }
    }

    return true;
}

int Polygon::shrink(Edge* edge) {
    auto it = find(edges_.begin(), edges_.end(), edge);
    if (it != edges_.end()) {
        edges_.erase(it);
    } else {
        printf("edge %ld not found in polygon %ld\n", edge->id_, id_);
        exit(1);
    }
    if (edges_.size() == 3) {
        for (auto e : edges_) {
            e->triangle_count_ = e->triangle_count_ + 1;
        }
    }
    if (edges_.size() == 2) {
        for (auto e : edges_) {
            e->triangle_count_ = e->triangle_count_ - 1;
        }
    }

    return 0;
}

int Polygon::expand(Edge* edge) {
    edges_.push_back(edge);
    if (edges_.size() == 3) {
        for (auto e : edges_) {
            e->triangle_count_ = e->triangle_count_ + 1;
        }
    }
    if (edges_.size() == 4) {
        for (auto e : edges_) {
            e->triangle_count_ = e->triangle_count_ - 1;
        }
    }

    return 0;
}

int Polygon::logEdges(std::string name) {
    printf("%s %ld\n",name.c_str(), id_);
    printf("%ld",edges_.size());
    for (auto edge : edges_) {
        printf(" %ld",edge->id_);
    }
    printf("\n");

    return 0;
}

int Polygon::SpringGon(double curtime) {
    int doneyet=0;
    if((type_==38||type_==94)){
    //if((type_==10000)){
        for (int i = 0; i < edges_.size(); i++) {
            double tempx = 0.0;
            double tempy = 0.0;
            double tempz = 0.0;
            double templ = 0.0;
            double tempforce = 0.0;
            int triggerdelam = 0;
            double tempSpringArea = 1.0;

            double ks = 1.0;
            tempx = edges_[i]->vertices_[0]->position_[0]-edges_[i]->vertices_[1]->position_[0];
            tempy = edges_[i]->vertices_[0]->position_[1]-edges_[i]->vertices_[1]->position_[1];
            tempz = edges_[i]->vertices_[0]->position_[2]-edges_[i]->vertices_[1]->position_[2];

            if(tempx > run_->Lx_/2.0){tempx += -run_->Lx_;}
            if(tempx < -run_->Lx_/2.0){tempx += run_->Lx_;}
            if(tempy > run_->Ly_/2.0){tempy += -run_->Ly_;}
            if(tempy < -run_->Ly_/2.0){tempy += run_->Ly_;}
            if(tempz > run_->Lz_/2.0){tempz += -run_->Lz_;}
            if(tempz < -run_->Lz_/2.0){tempz += run_->Lz_;}


            templ = pow(pow(tempx,2)+pow(tempy,2)+pow(tempz,2),0.5);

            //if(curtime > 151 && curtime < 151 + run_->dt_){std::cout << templ << std::endl;}
            //if(curtime > 200 && curtime < 200 + run_->dt_){std::cout << templ << std::endl;}

            for (int k = 0; k < 2; k++) {
                for (int m = 0; m < edges_[i]->vertices_[k]->cells_.size(); m++) {
                    if(edges_[i]->vertices_[k]->cells_[m]->type_==63){triggerdelam=1;
                    tempSpringArea=edges_[i]->vertices_[k]->cells_[m]->springArea_;
                //cout << tempSpringArea << endl;
                }
                }
            }

            tempforce = run_->mu_ * ks*(templ-run_->l0s_);
            if(triggerdelam==1){
                //if(templ-run_->l0s_ > 0 && templ < 0.3){
                if(templ-run_->l0s_ > 0 && tempforce < 0.5*(tempSpringArea)){
                edges_[i]->vertices_[0]->velocity_[0] += run_->mu_ * ks*(templ-run_->l0s_) * (-tempx/templ);
                edges_[i]->vertices_[1]->velocity_[0] += run_->mu_ * ks*(templ-run_->l0s_) * (tempx/templ);
                edges_[i]->vertices_[0]->velocity_[1] += run_->mu_ * ks*(templ-run_->l0s_) * (-tempy/templ);
                edges_[i]->vertices_[1]->velocity_[1] += run_->mu_ * ks*(templ-run_->l0s_) * (tempy/templ);
                edges_[i]->vertices_[0]->velocity_[2] += run_->mu_ * ks*(templ-run_->l0s_) * (-tempz/templ);
                edges_[i]->vertices_[1]->velocity_[2] += run_->mu_ * ks*(templ-run_->l0s_) * (tempz/templ);

                edges_[i]->vertices_[0]->springforces_[0] += run_->mu_ * ks*(templ-run_->l0s_) * (-tempx/templ);
                edges_[i]->vertices_[1]->springforces_[0] += run_->mu_ * ks*(templ-run_->l0s_) * (tempx/templ);
                edges_[i]->vertices_[0]->springforces_[1] += run_->mu_ * ks*(templ-run_->l0s_) * (-tempy/templ);
                edges_[i]->vertices_[1]->springforces_[1] += run_->mu_ * ks*(templ-run_->l0s_) * (tempy/templ);
                edges_[i]->vertices_[0]->springforces_[2] += run_->mu_ * ks*(templ-run_->l0s_) * (-tempz/templ);
                edges_[i]->vertices_[1]->springforces_[2] += run_->mu_ * ks*(templ-run_->l0s_) * (tempz/templ);


                }
            }

            else{
            edges_[i]->vertices_[0]->velocity_[0] += run_->mu_ * ks*(templ-run_->l0s_) * (-tempx/templ);
            edges_[i]->vertices_[1]->velocity_[0] += run_->mu_ * ks*(templ-run_->l0s_) * (tempx/templ);
            edges_[i]->vertices_[0]->velocity_[1] += run_->mu_ * ks*(templ-run_->l0s_) * (-tempy/templ);
            edges_[i]->vertices_[1]->velocity_[1] += run_->mu_ * ks*(templ-run_->l0s_) * (tempy/templ);
            edges_[i]->vertices_[0]->velocity_[2] += run_->mu_ * ks*(templ-run_->l0s_) * (-tempz/templ);
            edges_[i]->vertices_[1]->velocity_[2] += run_->mu_ * ks*(templ-run_->l0s_) * (tempz/templ);    

            edges_[i]->vertices_[0]->springforces_[0] += run_->mu_ * ks*(templ-run_->l0s_) * (-tempx/templ);
            edges_[i]->vertices_[1]->springforces_[0] += run_->mu_ * ks*(templ-run_->l0s_) * (tempx/templ);
            edges_[i]->vertices_[0]->springforces_[1] += run_->mu_ * ks*(templ-run_->l0s_) * (-tempy/templ);
            edges_[i]->vertices_[1]->springforces_[1] += run_->mu_ * ks*(templ-run_->l0s_) * (tempy/templ);
            edges_[i]->vertices_[0]->springforces_[2] += run_->mu_ * ks*(templ-run_->l0s_) * (-tempz/templ);
            edges_[i]->vertices_[1]->springforces_[2] += run_->mu_ * ks*(templ-run_->l0s_) * (tempz/templ);

            }
        }

    }


    return 0;
}

int Polygon::UpDateTension() {

    double curForceX = 0;
    double curForceY = 0;
    double curForceZ = 0;
    for (int i = 0; i < vertices_.size(); i++) { 
        curForceX += (vertices_[i]->springforces_[0]);
        curForceY += (vertices_[i]->springforces_[1]);
        curForceZ += (vertices_[i]->springforces_[2]);

        //curForceX += (vertices_[i]->velocity_[0]);
        //curForceY += (vertices_[i]->velocity_[1]);
        //curForceZ += (vertices_[i]->velocity_[2]);

        // if(type_==94){
        //     cout << i << endl;
        //     cout << vertices_[i]->springforces_[2] << endl;
        //     cout << curForceZ << endl;
        //}
    }

    springtension_[0] = curForceX;
    springtension_[1] = curForceY;
    springtension_[2] = curForceZ;

    return 0;
}

int Polygon::ResetTension() {

    springtension_[0] = 0;
    springtension_[1] = 0;
    springtension_[2] = 0;

    return 0;
}