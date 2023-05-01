// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <unordered_map>
#include <random>
#include <set>
#include "Run.h"

#include <experimental/filesystem>
#include <exception>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
//#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkIdList.h>
#include <vtkNamedColors.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyhedron.h>
#include <vtkPolygon.h>
//#include <vtkProperty.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersion.h>
#include <vtkXMLUnstructuredGridWriter.h>

using namespace std;

Run::Run() {
//    dt_ = 0.001;
//    dtr_ = 10*dt_;
//    dump_period_ = 10000*dt_; //
     log_period_ = 5000;
     t_start_ = 0.;
     t_end_ = 50000.;
    mu_ = 1.0;
    kB_ = 1.0;
//    temperature_ = 1.0e-5;
    //xxx
    Lx_ = 12.;
    Ly_ = 12.;
    Lz_ = 12.;
    NCell_ = 1728;
}
int     Run::InitializeCellUniform() {
    updateCellVertices();
    for (auto cell : cells_) {
        cell->type_ = 1;
        cell->s0_=5.6;
        cell->v0_=1.0;
    }
    
    return 0;
}

int     Run::InitializeCellDirectors() {
    std::default_random_engine generator(std::random_device{}());
    std::uniform_real_distribution<> disttheta(-1, 1);
    std::uniform_real_distribution<> distphi(0, 2*M_PI);
    double temptheta;
    for (auto cell : cells_) {
        temptheta = disttheta(generator);
        cell->cellDirectors_[0]=acos(temptheta);
        cell->cellDirectors_[1]=distphi(generator);
    }
    
    return 0;
}

int     Run::InitializeCellDelam() {
    int check = 0;

    int delam = 1;
    int setup = 0;

    if(delam==1){
    for (auto cell : cells_) {
    	int check1=0;
    	int check2=0;
    	if(check==0)
    	{
    		for (auto polygon : cell->polygons_) {
    			//2 and 4/16 for some reason
    			//if(polygon->type_==7 || polygon->type_== 38){check1=1;}
          if(polygon->type_==7){check1=1;}
    			if(polygon->type_==0){check2=1;}
          //if(polygon->type_==22){check2=1;}
    		}
    		if(check1==1 && check2==1){
		        cell->type_ = 63;
		        cell->s0_=s07_;
		        cell->v0_=v07_;
		    	check=1;
		    	cell->color_ = 3;
		    }
        }
        if(cell->type_==3){
        cell->type_ = 15;
        cell->s0_=s05_;
        cell->v0_=v05_;
        }        
        if(cell->type_==1){
        cell->type_ = 7;
        cell->s0_=s04_;
        cell->v0_=v04_;
        }    
    }}

    if(setup==1){
        //temperature_ = 0.01;
        for (auto cell : cells_) {
            if(cell->type_==1){
            //cell->type_ = 7;
            cell->s0_=s04_;
            cell->v0_=v04_;
            }
            if(cell->type_==3){
            //cell->type_ = 15;
            cell->s0_=s05_;
            cell->v0_=v05_;
            }
        }
    }
    
    return 0;
}

int     Run::InitializePlacode() {
    int check = 0;

    for (auto cell : cells_) {
    	if(cell->type_==1){
    		if(cell->cellcenter_[0] > Lx_/2-2 && cell->cellcenter_[0] < Lx_/2+2 && (cell->cellcenter_[1]<2 || cell->cellcenter_[1]>Ly_-2)){
	        cell->type_ = 63;
          cell->color_ = 3;
	        cell->s0_=s07_;
	        cell->v0_=v07_;
	    	check+=1;
        	}
    	}
        if(cell->type_==3){
        cell->type_ = 15;
        cell->s0_=s05_;
        cell->v0_=v05_;
        }        
        if(cell->type_==1){
        cell->type_ = 7;
        cell->s0_=s04_;
        cell->v0_=v04_;
        } 
    }
    
    return 0;
}

int     Run::PlacodeCenter() {
  double countplac = 0;
  double tempx = 0;
  double tempy = 0;
  double tempthx  = 0;
  double tempphix  = 0;
  double tempthy = 0;
  double tempphiy  = 0;
  for (auto cell : cells_){
    if(cell->type_ == 63){
      countplac += 1.0;
      tempx = 2*M_PI*cell->cellcenter_[0]/Lx_;
      tempy = 2*M_PI*cell->cellcenter_[1]/Ly_;

      tempthx += cos(tempx);
      tempphix += sin(tempx);
      tempthy += cos(tempy);
      tempphiy += sin(tempy);

    }
  }
  tempthx = tempthx/countplac;
  tempthy = tempthy/countplac;
  tempphix = tempphix/countplac;
  tempphiy = tempphiy/countplac;

  placodex_ = Lx_*(atan2(-tempphix,-tempthx)+M_PI)/(2*M_PI);
  placodey_ = Ly_*(atan2(-tempphiy,-tempthy)+M_PI)/(2*M_PI);
  placodesize_ = countplac;
  cout << placodex_ << endl;
  cout << placodey_ << endl;
  
  return 0;
}


int     Run::InitalizeErin() {
    int randomcell = rand() % NCell_;
    cout << randomcell << endl;
    cells_[randomcell]->type_ = 1;
    cells_[randomcell]->s0_=s02_;
    cells_[randomcell]->v0_=v02_;    
    
    return 0;
}

int     Run::InitalizeSolidifcation() {
    int check = 0;

    int setup = 1;

    if(setup==1){
        for (auto cell : cells_) {
          if(cell->type_==0){
          cell->type_ = 31;         
          }   
        }
    }
    
    return 0;
}


int     Run::InitializeCellType() {
    //type1count = NCell_/2;
    //type2count = NCell_-type1count;
    type1count = NCell_/2;
    type2count = NCell_/2;
    int randompos = 0;
    int topbot =1;
    int stratified = 0; //ECM network
    int nbasment = 0;
    int nbasal = 0;
    int nsupra = 0;

    updateCellVertices();
    if(randompos==1 || stratified==1){
    for (auto cell : cells_) {
        if(randompos==1)
        {
            if(rand() % (type1count+type2count) < type1count){
                cell->type_ = 0;
                cell->color_ = 0;
                type1count -=1;
            }
            else{
                cell->type_ = 7;
                cell->color_ = 1;
                type2count -=1;           
            }
        }
        if(stratified==1){
        double layer = 0;
        //int type1countmax = Lx_*Ly_-round(Lx_*Ly_*0.05);
        int type1countmax = Lx_*Ly_+round(Lx_*Ly_*0.05);
        type1count = type1countmax;
        double layermin = Lz_/2;

        while(type1count>0)
        {
            layer += 0.01;
            //cout << layer << endl;
            type1count = type1countmax;
            nsupra = 0;
            nbasment = 0;
            nbasal = 0;


            for (auto cell : cells_) 
            {
                if(type1count>0)
                {   
                   //cout << layer*Lz_/10.0 << endl;
                   if(cell->cellcenter_[2]<layermin + layer*(Lz_-layermin)/10.0 && cell->cellcenter_[2]>layermin)
                    {
                        cell->type_ = 7;
                        cell->color_ = 1;
                        type1count -=1;
                        nbasal +=1;
                    }
                    else if(cell->cellcenter_[2]<layermin){
                        cell->type_ = 0;
                        cell->color_ = 0;
                        nbasment +=1;
                    }
                    else
                    {
                        cell->type_ = 15;
                        cell->color_ = 2;     
                        nsupra +=1;      
                    }                 
                }
            }
        }   
    //cout << nbasment << endl;
    //cout << nbasal << endl;
    //cout << nsupra << endl;
    }


    if(cell->type_ == 0){
      cell->s0_=s01_;
      //cell->s0_=6.65;
      cell->v0_=v01_;
      }
    if(cell->type_ == 7){
      cell->s0_=s04_;
      //cell->s0_=6.65;
      cell->v0_=v04_;
      }
    if(cell->type_ == 15){
      cell->s0_=s05_;
      //cell->s0_=6.65;
      cell->v0_=v05_;
      }

    }}

    if(topbot==1)
    {
        double layer = 0;
        type1count = NCell_/2;
        while(type1count>0)
        {
            layer += 0.01;
            //cout << layer << endl;
            type1count = NCell_/2;
            for (auto cell : cells_) 
            {
                if(type1count>0)
                {   
                   //cout << layer*Lz_/10.0 << endl;
                   if(cell->cellcenter_[2]<layer*Lz_/10.0)
                    {
                        cell->type_ = 0;
                        cell->color_ = 0;
                        type1count -=1;
                    }
                    else
                    {
                        cell->type_ = 7;     
                        cell->color_ = 1;      
                    }                 
                }
            }
        }         
    for (auto cell : cells_) 
        {   
        if(cell->type_ == 0){
          cell->s0_=s01_;
          cell->v0_=v01_;
          }
        if(cell->type_ == 1){
          cell->s0_=s02_;
          cell->v0_=v02_;
          }
        if(cell->type_ == 3){
          cell->s0_=s03_;
          cell->v0_=v03_;
          }
    if(cell->type_ == 7){
      cell->s0_=s04_;
      //cell->s0_=6.65;
      cell->v0_=v04_;
      }        
        }
    }

    
    if(stratified==1){
      cout << "Number of Basement Cells: " << nbasment << endl;
      cout << "Number of Basal Cells: " << nbasal << endl;
      cout << "Number of Suprabasal Cells: " << nsupra << endl;
    }
    else{
    cout << type1count << endl;
  }


    return 0;
}


int Run::start() {
    count_reconnect_ = 0;
    count_dump_ = 0;
    count_log_ = 0;
    simulation_time_ = t_start_;
    double t_roundError = 0.01*dt_;
    auto start = chrono::steady_clock::now();

    printf("\nSimulation Start ...\n");
    printf("Real time elapsed: Rte\n");
    printf("Time        ");
    printf("Rte         ");
    printf("Volume      ");
    printf("I->H        ");
    printf("H->I        ");
    printf("E_volume    ");
    printf("E_interface ");
    //printf("Type of cell 0");
    printf("Energy ");
    printf("Average Force ");
    printf("Max Force      \n");


      // ------------------------------------------
      // Paraview output stuff
      // ------------------------------------------

      // choose whether to output to Paraview
      bool outputVtkVis = true;

      // check if paraview output directory exists, if not, create it
      if( outputVtkVis ) {
        if( !std::experimental::filesystem::exists( "paraview" ) ) {
          std::experimental::filesystem::create_directory("paraview");
        }
        else {
          // check if is indeed a directory (it could be a file)
          if( !std::experimental::filesystem::is_directory("paraview")) {
            std::cerr << "\"paraview\" is a file, not a directory" << std::endl;
            std::cerr << "Remove \"paraview\" file, and re-run the program." << std::endl;
            throw std::exception();
          }
        }
      }

      // header of timeseries file
      ofstream vtkTimeseries ("timeseries.pvd");
      if( outputVtkVis ) {
        vtkTimeseries << "<?xml version=\"1.0\"?>" << endl;
        vtkTimeseries << "<VTKFile type=\"Collection\" version=\"0.1\"" << endl;
        vtkTimeseries << "         byte_order=\"LittleEndian\"" << endl;
        vtkTimeseries << "         compressor=\"vtkZLibDataCompressor\">" << endl;
        vtkTimeseries << "  <Collection>" << endl;
      }


    while (simulation_time_ < t_end_ + t_roundError) {
    	
        //Old
    	//if(simulation_time_ > dt_ && simulation_time_ < 5*dt_){InitializeCellType();}
    	
        // After letting system iniatilize apply actual cell types:
    	//if(simulation_time_ > log_period_ + dt_ && simulation_time_ < log_period_ + 2*dt_){InitializeCellType();}
    	//if(simulation_time_ > 2*log_period_ + dt_ && simulation_time_ < 2*log_period_ + 2*dt_){InitializeCellDelam();}
    	

    	if(simulation_time_ > 25 && simulation_time_ < 25 + dt_){InitializeCellType();
            //vertices_[0]->updateSP(temperature_);
        }
       if(simulation_time_ > 75 && simulation_time_ < 75 + dt_){temperature_ = 0.0;}
       

        //if(simulation_time_ > 75 + dt_ && simulation_time_ < 75 + 2*dt_){InitalizeSolidifcation();}
        if(simulation_time_ > 200 && simulation_time_ < 200 + dt_){
          //temperature_=0.0;
          InitializeCellDelam();}
        

    	  //if(simulation_time_ > 100 + dt_ && simulation_time_ < 100 + 2*dt_){InitializePlacode();
          //placodeon_= 1;
      	//}
        //if(placodeon_==1){PlacodeCenter();}
        
        // update geometry information
        updateGeoinfo();
        // update volumeForces
        volume_->updateForces();
        // update interfaceForces
        interface_->updateForces();
        // update velocities
        updateVerticesVelocity();
        updatePolygonVertices();

        //Update cell interface type and cell directors
        for (auto cell : cells_) {
            for (auto polygon : cell->polygons_) {
                polygon->type_ = 0;
            }
        }

        for (auto cell : cells_) {            
            for (auto polygon : cell->polygons_) {
                polygon->type_ += cell->type_;
            }
        }        


        // log to screen
        if (simulation_time_ - t_start_ + t_roundError  > count_log_ * log_period_ || simulation_time_>200) {
            double demixing = 0.;
            double countcells = 0.;

            double totalforce = 0.0;
            double maxforce = 0.0;
            for (long int i = 0; i < vertices_.size(); i++) {
              vertices_[i]->magForce();
              totalforce += vertices_[i]->netforce_;
              if(vertices_[i]->netforce_>maxforce){
                maxforce = vertices_[i]->netforce_;
              }
            }
            totalforce = totalforce/vertices_.size();

            volume_->updateEnergy();
            interface_->updateEnergy();
            printf("%-12.2f%-12.3f%-12.3f%-12ld%-12ld%-12.6f%-12.6f%-12.6f%-12.6f%-12.6f\n", simulation_time_,
                   (chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count())/1.0e6,
                   volume_->totalVolume_,
                   reconnection_->count_IH_,
                   reconnection_->count_HI_,
                   volume_->energy_,
                   interface_->energy_,
                   volume_->energy_+interface_->energy_,
                   totalforce,
                   maxforce);
            start = chrono::steady_clock::now();
            reconnection_->count_IH_ = 0;
            reconnection_->count_HI_ = 0;
            count_log_++;

            //Experimental Paraview addition
            if(outputVtkVis ) {
               // Paraview file name
             char filename[256]; 
             sprintf(filename,"paraview/cells_t0_%06f.vtu", simulation_time_);

             //int nCells = t.cells().size();
             //std::cout << "Number of cells " << NCell_ << std::endl;
             
             // unstructured grid and points (i.e. vertices)
             vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
             vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

             updatePolygonVertices();
             long int Npolygons = polygons_.size();

             // time stamp
             vtkSmartPointer<vtkDoubleArray> vtkTimeArray = vtkSmartPointer<vtkDoubleArray>::New();
             vtkTimeArray->SetNumberOfComponents(1);
             vtkTimeArray->SetNumberOfTuples(1);
             vtkTimeArray->SetName("TIME");            
 
             // polygon ID
             vtkSmartPointer<vtkIntArray> polyID = vtkSmartPointer<vtkIntArray>::New();
             polyID->SetNumberOfComponents(1);
             polyID->SetNumberOfTuples(polygons_.size());
             polyID->SetName("polyID");

             // polygon type
             vtkSmartPointer<vtkIntArray> polyType = vtkSmartPointer<vtkIntArray>::New();
             polyType->SetNumberOfComponents(1);
             polyType->SetNumberOfTuples(polygons_.size());
             polyType->SetName("polyType");

             // cellID
             vtkSmartPointer<vtkIntArray> cellID = vtkSmartPointer<vtkIntArray>::New();
             cellID->SetNumberOfComponents(1);
             cellID->SetNumberOfTuples(polygons_.size());
             cellID->SetName("cellID");

             // cellType
             vtkSmartPointer<vtkIntArray> cellType = vtkSmartPointer<vtkIntArray>::New();
             cellType->SetNumberOfComponents(1);
             cellType->SetNumberOfTuples(polygons_.size());
             cellType->SetName("cellType");

             // cell center
             vtkSmartPointer<vtkDoubleArray> cellCenter = vtkSmartPointer<vtkDoubleArray>::New();
             cellCenter->SetNumberOfComponents(3);
             cellCenter->SetNumberOfTuples(polygons_.size());
             cellCenter->SetName("cellCenter");

             // Polygon neighbors
             vtkSmartPointer<vtkDoubleArray> PolyNeighs = vtkSmartPointer<vtkDoubleArray>::New();
             PolyNeighs->SetNumberOfComponents(2);
             PolyNeighs->SetNumberOfTuples(polygons_.size());
             PolyNeighs->SetName("PolyNeighs");

             // cell volume
             vtkSmartPointer<vtkDoubleArray> cellVol = vtkSmartPointer<vtkDoubleArray>::New();
             cellVol->SetNumberOfComponents(1);
             cellVol->SetNumberOfTuples(polygons_.size());
             cellVol->SetName("cellVol");

             // cell SA
             vtkSmartPointer<vtkDoubleArray> cellArea = vtkSmartPointer<vtkDoubleArray>::New();
             cellArea->SetNumberOfComponents(1);
             cellArea->SetNumberOfTuples(polygons_.size());
             cellArea->SetName("cellArea");

             // Polygon Area
             vtkSmartPointer<vtkDoubleArray> polyArea = vtkSmartPointer<vtkDoubleArray>::New();
             polyArea->SetNumberOfComponents(1);
             polyArea->SetNumberOfTuples(polygons_.size());
             polyArea->SetName("polyArea");

             // Polygon Force
             vtkSmartPointer<vtkDoubleArray> polyForce = vtkSmartPointer<vtkDoubleArray>::New();
             polyForce->SetNumberOfComponents(3);
             polyForce->SetNumberOfTuples(polygons_.size());
             polyForce->SetName("polyForce");

	         // //cell registration general
	         // vtkSmartPointer<vtkDoubleArray> cellRegister12 = vtkSmartPointer<vtkDoubleArray>::New();
	         // cellRegister12->SetNumberOfComponents(1);
	         // cellRegister12->SetNumberOfTuples(polygons_.size());
	         // cellRegister12->SetName("cellRegister12");

          //    // cell ellipsoid principal a-axis (major)
          //    vtkSmartPointer<vtkDoubleArray> cellSP = vtkSmartPointer<vtkDoubleArray>::New();
          //    cellSP->SetNumberOfComponents(3);
          //    cellSP->SetNumberOfTuples(polygons_.size());
          //    cellSP->SetName("cellSP");

          //    // cell ellipsoid principal a-axis (major)
          //    vtkSmartPointer<vtkDoubleArray> cellSPa = vtkSmartPointer<vtkDoubleArray>::New();
          //    cellSPa->SetNumberOfComponents(2);
          //    cellSPa->SetNumberOfTuples(polygons_.size());
          //    cellSPa->SetName("cellSPa");

             // cell ellipsoid principal a-axis (major)
             vtkSmartPointer<vtkDoubleArray> cellEllipsoidAaxis = vtkSmartPointer<vtkDoubleArray>::New();
             cellEllipsoidAaxis->SetNumberOfComponents(3);
             cellEllipsoidAaxis->SetNumberOfTuples(polygons_.size());
             cellEllipsoidAaxis->SetName("cellEllipsoidAaxis");

             // cell ellipsoid principal b-axis 
             vtkSmartPointer<vtkDoubleArray> cellEllipsoidBaxis = vtkSmartPointer<vtkDoubleArray>::New();
             cellEllipsoidBaxis->SetNumberOfComponents(3);
             cellEllipsoidBaxis->SetNumberOfTuples(polygons_.size());
             cellEllipsoidBaxis->SetName("cellEllipsoidBaxis");

             // cell ellipsoid principal c-axis (minor)
             vtkSmartPointer<vtkDoubleArray> cellEllipsoidCaxis = vtkSmartPointer<vtkDoubleArray>::New();
             cellEllipsoidCaxis->SetNumberOfComponents(3);
             cellEllipsoidCaxis->SetNumberOfTuples(polygons_.size());
             cellEllipsoidCaxis->SetName("cellEllipsoidCaxis");

             // cell ellipsoid principal a-axis orientation
	         vtkSmartPointer<vtkDoubleArray> cellorient = vtkSmartPointer<vtkDoubleArray>::New();
	         cellorient->SetNumberOfComponents(1);
	         cellorient->SetNumberOfTuples(polygons_.size());
	         cellorient->SetName("cellOrientation");

             // cell anisotropy
	         vtkSmartPointer<vtkDoubleArray> cellani = vtkSmartPointer<vtkDoubleArray>::New();
	         cellani->SetNumberOfComponents(1);
	         cellani->SetNumberOfTuples(polygons_.size());
	         cellani->SetName("cellAnisotropy");


	         // cell ellipsoid principal radii
	         vtkSmartPointer<vtkDoubleArray> cellEllipsoidRadii = vtkSmartPointer<vtkDoubleArray>::New();
	         cellEllipsoidRadii->SetNumberOfComponents(3);
	         cellEllipsoidRadii->SetNumberOfTuples(polygons_.size());
	         cellEllipsoidRadii->SetName("cellEllipsoidRadii");

            //cell registration z
            vtkSmartPointer<vtkDoubleArray> cellcolour = vtkSmartPointer<vtkDoubleArray>::New();
            cellcolour->SetNumberOfComponents(1);
            cellcolour->SetNumberOfTuples(polygons_.size());
            cellcolour->SetName("cellcolour");


            int polyCounter = 0;
            int polyCounter2 = 0;
            int polyCounter3 = 0;

            int numberTotalVerticesWithDups = 0;

            int pointtotal = 0;
            int cellcounter = 0;


            for (long int i = 0; i < polygons_.size(); i++) {
            	polygons_[i]->twincell_[0] = -1;
            	polygons_[i]->twincell_[1] = -1; 
            	polygons_[i]->twincell_[2] = -1; 
            	polygons_[i]->twincell_[3] = -1; 
            }
            cout << "Total poly: " << polygons_.size() << endl;

            for (auto cell : cells_) {
	            int xplus = 0;
	            int xminus = 0;
	            int yplus = 0;
	            int yminus = 0;
	            int zplus = 0;
	            int zminus =0;


            	for (long int i = 0; i < cell->polygons_.size(); i++) { 



	            	for (int j = 0; j < cell->polygons_[i]->vertices_.size(); j++) { 
		                

		                int numberOfVerticesOfFace = cell->polygons_[i]->vertices_.size();

		                double polyverts[numberOfVerticesOfFace][3];
	                    polyverts[j][0]=cell->polygons_[i]->vertices_[j]->position_[0];
	                    polyverts[j][1]=cell->polygons_[i]->vertices_[j]->position_[1];
	                    polyverts[j][2]=cell->polygons_[i]->vertices_[j]->position_[2];

	                    if(polyverts[j][0]<1){xminus=1;}
	                    if(polyverts[j][0]>Lx_-1){xplus=1;}
	                    if(polyverts[j][1]<1){yminus=1;}
	                    if(polyverts[j][1]>Ly_-1){yplus=1;}
	                    if(polyverts[j][2]<1){zminus=1;}
	                    if(polyverts[j][2]>Lz_-1){zplus=1;}
            		}
            	}

              double tempVec[3];
              double totVec[3]{0,0,0};
              double totArea = 0;
              double totalvert[3];
              double tempArea = 0;
              double cell0neighs[16]{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
              int cell0neighcount = 0;
              int counttemp = 0;

		    	for (long int i = 0; i < cell->polygons_.size(); i++) {  
                polyCounter++;


	                vtkSmartPointer<vtkIdList> vtkFaces = vtkSmartPointer<vtkIdList>::New(); // vtk faces  
	                //vtkFaces->InsertNextId(polygons_[i]->vertices_.size());
	                
	                int numberOfVerticesOfFace = cell->polygons_[i]->vertices_.size();
	                numberTotalVerticesWithDups += cell->polygons_[i]->vertices_.size();

	                double polyverts[numberOfVerticesOfFace][3];

                  tempVec[0]=0;
                  tempVec[1]=0;
                  tempVec[2]=0;
                  counttemp = 0;
                  tempArea = cell->polygons_[i]->area_;
                  totArea += cell->polygons_[i]->area_;
	                
	                //cout << "Number of Verticies of Polygon: " << polygons_[i]->vertices_.size() << endl;
	                for (int j = 0; j < cell->polygons_[i]->vertices_.size(); j++) {

	                    polyverts[j][0]=cell->polygons_[i]->vertices_[j]->position_[0];
	                    polyverts[j][1]=cell->polygons_[i]->vertices_[j]->position_[1];
	                    polyverts[j][2]=cell->polygons_[i]->vertices_[j]->position_[2];

	                    if(xminus==1 && xplus==1 && polyverts[j][0]<2){polyverts[j][0]+=Lx_;}
	                    if(yminus==1 && yplus==1 && polyverts[j][1]<2){polyverts[j][1]+=Ly_;}
	                    if(zminus==1 && zplus==1 && polyverts[j][2]<2){polyverts[j][2]+=Lz_;}

                      tempVec[0] += polyverts[j][0];
                      tempVec[1] += polyverts[j][1];
                      tempVec[2] += polyverts[j][2];
                      counttemp +=1;

	                    points->InsertNextPoint(polyverts[j][0],polyverts[j][1],polyverts[j][2]); 
	                    //cout << "Vertex: " << pointtotal << " Estimated: " << numberTotalVerticesWithDups-polygons_[i]->vertices_.size()+j << " X " << polyverts[j][0] << " Y " << polyverts[j][1] << " Z " << polyverts[j][2] << endl;
	                    pointtotal++;
	                    }             
	                //cout << "Done" << endl;
                totVec[0] += tempVec[0]*tempArea/counttemp;
                totVec[1] += tempVec[1]*tempArea/counttemp;
                totVec[2] += tempVec[2]*tempArea/counttemp;

                //vtkFaces->InsertNextId(polygons_[i]->vertices_.size());
                for( int jj=0 ; jj < cell->polygons_[i]->vertices_.size() ; ++jj ) { 
                  vtkFaces->InsertNextId(numberTotalVerticesWithDups-cell->polygons_[i]->vertices_.size()+jj);
                    }    

                double tempnorm[3];
                
                tempnorm[0]=0;
                tempnorm[1]=0;
                tempnorm[2]=0;
                int currpos;
                int nextpos;


                for (int j = 0; j < cell->polygons_[i]->vertices_.size(); j++) {
                	currpos = j;
                	nextpos = j+1;

                  if(j+1==cell->polygons_[i]->vertices_.size()){nextpos = 0;}



                	tempnorm[0] += (polyverts[currpos][1]-polyverts[nextpos][1])*(polyverts[currpos][2]+polyverts[nextpos][2]);
                	tempnorm[1] += (polyverts[currpos][2]-polyverts[nextpos][2])*(polyverts[currpos][0]+polyverts[nextpos][0]);
                	tempnorm[2] += (polyverts[currpos][0]-polyverts[nextpos][0])*(polyverts[currpos][1]+polyverts[nextpos][1]);

                }

                double tempmagvec;
                tempmagvec = pow(pow(tempnorm[0],2)+pow(tempnorm[1],2)+pow(tempnorm[2],2),0.5);
                cell->polygons_[i]->normvector_[0] = tempnorm[0]/tempmagvec;
                cell->polygons_[i]->normvector_[1] = tempnorm[1]/tempmagvec;
                cell->polygons_[i]->normvector_[2] = tempnorm[2]/tempmagvec;


	            // add cell to unstructure grid
	            uGrid->InsertNextCell(VTK_POLYGON,vtkFaces);
	            polyID->InsertValue(polyCounter-1, polyCounter);
	            polyType->InsertValue(polyCounter-1, cell->polygons_[i]->type_);

              polyArea->InsertValue(polyCounter-1, cell->polygons_[i]->area_);
              polyForce->InsertTuple(polyCounter-1, cell->polygons_[i]->springtension_);



	            uGrid->GetCellData()->AddArray(polyType);
	            uGrid->GetCellData()->AddArray(polyID);
              uGrid->GetCellData()->AddArray(polyArea);
              uGrid->GetCellData()->AddArray(polyForce);
              }
              

              totVec[0] = totVec[0]/totArea;
              totVec[1] = totVec[1]/totArea;
              totVec[2] = totVec[2]/totArea;

              //cout << "Current Npoly: " << cell->polygons_.size() << endl;
              for (long int i = 0; i < cell->polygons_.size(); i++) { 
              polyCounter2++; 
              int curcellID = cell->id_;
              int curpolyID = cell->polygons_[i]->id_;

              if(cell->polygons_[i]->twincell_[0] ==-1){cell->polygons_[i]->twincell_[0] = cell->id_;
                  cell->polygons_[i]->twincell_[2] = curpolyID;}
              else if(cell->polygons_[i]->twincell_[2] == curpolyID){cell->polygons_[i]->twincell_[1] = cell->id_;
                  cell->polygons_[i]->twincell_[3] = curpolyID;}
            else{cout << "Triple Cell?" << endl;}

              cell->polygons_[i]->cellcenter_[0] = totVec[0];
              cell->polygons_[i]->cellcenter_[1] = totVec[1];
              cell->polygons_[i]->cellcenter_[2] = totVec[2];

              cells_[cell->id_]->cellcenter_[0]  = totVec[0];
              cells_[cell->id_]->cellcenter_[1]  = totVec[1];
              cells_[cell->id_]->cellcenter_[2]  = totVec[2];             

              double tempdirect[3];
              tempdirect[0] = cos(cell->cellDirectors_[1])*sin(cell->cellDirectors_[0]);
              tempdirect[1] = sin(cell->cellDirectors_[1])*sin(cell->cellDirectors_[0]);
              tempdirect[2] = cos(cell->cellDirectors_[0]);
              double tempangles[2];
              tempangles[0] = cell->cellDirectors_[0];
              tempangles[1] = cell->cellDirectors_[1];

              //cellSP->InsertTuple(polyCounter2-1, tempdirect);
              //cellSPa->InsertTuple(polyCounter2-1, tempangles);
              cellCenter->InsertTuple(polyCounter2-1, totVec); 
              cellID->InsertValue(polyCounter2-1, curcellID);
              cellArea->InsertValue(polyCounter2-1, totArea);
              cellVol->InsertValue(polyCounter2-1, cell->volume_);
              cellType->InsertValue(polyCounter2-1, cell->type_);
              cellcolour->InsertValue(polyCounter2-1, cell->color_);

              uGrid->GetCellData()->AddArray(cellID);
              //uGrid->GetCellData()->AddArray(cellSP);
              //uGrid->GetCellData()->AddArray(cellSPa);
              uGrid->GetCellData()->AddArray(cellCenter);
              uGrid->GetCellData()->AddArray(cellArea);
              uGrid->GetCellData()->AddArray(cellVol);
              uGrid->GetCellData()->AddArray(cellType);
              uGrid->GetCellData()->AddArray(cellcolour);
	            }
	       }

	     double distarray=0;
	     double xarray = 0;
	     double yarray = 0;
	     double zarray = 0;
       int check = 0;
       double vecmagave=0;
       int vecmagnum = 0;
       double xyave=0;
       int xyavenum = 0;
         for (auto cell : cells_) {
            EllipsoidByUnitPointMassPolyhedron fittedCell = fitEllipsoidByUnitPointMassPolyhedron(cell->id_);
            double tempAaxis[3];
            double tempBaxis[3];
            double tempCaxis[3];
            double tempradii[3];
            double temporient;
            double tempani;

            tempAaxis[0] = fittedCell.aAxis().x();
            tempAaxis[1] = fittedCell.aAxis().y();
            tempAaxis[2] = fittedCell.aAxis().z();

            tempBaxis[0] = fittedCell.bAxis().x();
            tempBaxis[1] = fittedCell.bAxis().y();
            tempBaxis[2] = fittedCell.bAxis().z();

            tempCaxis[0] = fittedCell.cAxis().x();
            tempCaxis[1] = fittedCell.cAxis().y();
            tempCaxis[2] = fittedCell.cAxis().z();

            temporient = acos(abs(tempAaxis[2]));

            tempradii[0] = fittedCell.a();
            tempradii[1] = fittedCell.b();
            tempradii[2] = fittedCell.c();
            
            tempani = (tempradii[2]/tempradii[0])*(tempradii[1]/tempradii[0]);


          for (long int i = 0; i < cell->polygons_.size(); i++) { 
                polyCounter3++;
                double tempreg=0;
                double tempregZ=0;
                double TotalCells=cells_.size();
                double tempx=0;
                double tempy=0;
                double tempz=0;
                double parallel=0;
                double perpendicular=0;
                double vecmag=0;
                double tempneighs[2];
                int curcellID = cell->id_;
                

                tempx = cells_[cell->polygons_[i]->twincell_[0]]->cellcenter_[0] - cells_[cell->polygons_[i]->twincell_[1]]->cellcenter_[0];
                tempy = cells_[cell->polygons_[i]->twincell_[0]]->cellcenter_[1] - cells_[cell->polygons_[i]->twincell_[1]]->cellcenter_[1];
                tempz = cells_[cell->polygons_[i]->twincell_[0]]->cellcenter_[2] - cells_[cell->polygons_[i]->twincell_[1]]->cellcenter_[2];

                if(tempx>Lx_/2){tempx = tempx - Lx_;}
                if(tempx<-Lx_/2){tempx = tempx + Lx_;}
                if(tempy>Ly_/2){tempy = tempy - Ly_;}
                if(tempy<-Ly_/2){tempy = tempy + Ly_;}
                if(tempz>Lz_/2){tempz = tempz - Lz_;}
                if(tempz<-Lz_/2){tempz = tempz + Lz_;}

                parallel = tempx*cell->polygons_[i]->normvector_[0]+tempy*cell->polygons_[i]->normvector_[1]+tempz*cell->polygons_[i]->normvector_[2];

                perpendicular  = pow(pow(tempx,2)+pow(tempy,2)+pow(tempz,2)-pow(parallel,2),0.5);
                vecmag = pow(pow(tempx,2)+pow(tempy,2)+pow(tempz,2),0.5);
                vecmagave+=vecmag;
                vecmagnum+=1;

                tempreg = perpendicular;

                tempx = abs(cells_[cell->polygons_[i]->twincell_[0]]->cellcenter_[0] - cells_[cell->polygons_[i]->twincell_[1]]->cellcenter_[0]);
                tempy = abs(cells_[cell->polygons_[i]->twincell_[0]]->cellcenter_[1] - cells_[cell->polygons_[i]->twincell_[1]]->cellcenter_[1]);


                //if(cell->polygons_[i]->area_ > 0.5){
                if(abs(tempx)>Lx_-2){
                  tempx=abs(tempx)-Lx_;
                }
                if(abs(tempy)>Ly_-2){
                  tempy=abs(tempy)-Ly_;
                }

                tempregZ = 1-pow(pow(tempx,2)+pow(tempy,2),0.5)/(1.53);
 
                xyave+=pow(pow(tempx,2)+pow(tempy,2),0.5);
                xyavenum+=1;

                distarray+=pow(pow(tempx,2)+pow(tempy,2)+pow(tempz,2),0.5);
                xarray+=abs(tempx);
                yarray+=abs(tempy);
                zarray+=abs(tempz);


                cellEllipsoidAaxis->InsertTuple(polyCounter3-1,tempAaxis);
                uGrid->GetCellData()->AddArray(cellEllipsoidAaxis);

                cellEllipsoidBaxis->InsertTuple(polyCounter3-1,tempBaxis);
                uGrid->GetCellData()->AddArray(cellEllipsoidBaxis);

                cellEllipsoidCaxis->InsertTuple(polyCounter3-1,tempCaxis);
                uGrid->GetCellData()->AddArray(cellEllipsoidCaxis);

                cellEllipsoidRadii->InsertTuple(polyCounter3-1, tempradii);   
                uGrid->GetCellData()->AddArray(cellEllipsoidRadii);

                cellorient->InsertValue(polyCounter3-1, temporient);   
                uGrid->GetCellData()->AddArray(cellorient);       

                cellani->InsertValue(polyCounter3-1, tempani);   
                uGrid->GetCellData()->AddArray(cellani);   

                //cellRegister12->InsertValue(polyCounter3-1, tempreg);
                //uGrid->GetCellData()->AddArray(cellRegister12);

                tempneighs[0] = cell->polygons_[i]->twincell_[0];
                tempneighs[1] = cell->polygons_[i]->twincell_[1];
                PolyNeighs->InsertTuple(polyCounter3-1, tempneighs); 
                uGrid->GetCellData()->AddArray(PolyNeighs);
              }
            }

             // add time stamp to vtk data set
             vtkTimeArray->InsertValue(0,simulation_time_);
             uGrid->GetFieldData()->AddArray(vtkTimeArray);

             uGrid->SetPoints(points);

             vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtkWriter =
                vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
             vtkWriter->SetInputData(uGrid);
             vtkWriter->SetFileName(filename);
             vtkWriter->SetDataModeToBinary();
             vtkWriter->Update();

              } // end outputVtkVis if

        }
        // dump
        if (simulation_time_ - t_start_ + t_roundError > count_dump_ * dump_period_) {
            if (simulation_time_ > (-0.01)*dt_) {
                dumpTopo();
                dumpCellCenter();
                dumpCellShapeIndex();
                dumpReconnection();
                dumpConfigurationVtk();

            }
//            dumpCellCenter();
//            dumpCellShapeIndex();
            count_dump_++;
        }

        // Euler dynamics
        updateVerticesPosition();

        // reconnect
        if (simulation_time_ - t_start_ + t_roundError > count_reconnect_ * dtr_) {
            reconnection_->start();
            count_reconnect_++;
        }


            for (auto cell : cells_) {
	            int xplus = 0;
	            int xminus = 0;
	            int yplus = 0;
	            int yminus = 0;
	            int zplus = 0;
	            int zminus =0;

            	for (long int i = 0; i < cell->polygons_.size(); i++) { 

	            	for (int j = 0; j < cell->polygons_[i]->vertices_.size(); j++) {                 

		                int numberOfVerticesOfFace = cell->polygons_[i]->vertices_.size();

		                double polyverts[numberOfVerticesOfFace][3];
	                    polyverts[j][0]=cell->polygons_[i]->vertices_[j]->position_[0];
	                    polyverts[j][1]=cell->polygons_[i]->vertices_[j]->position_[1];
	                    polyverts[j][2]=cell->polygons_[i]->vertices_[j]->position_[2];

	                    if(polyverts[j][0]<1){xminus=1;}
	                    if(polyverts[j][0]>Lx_-1){xplus=1;}
	                    if(polyverts[j][1]<1){yminus=1;}
	                    if(polyverts[j][1]>Ly_-1){yplus=1;}
	                    if(polyverts[j][2]<1){zminus=1;}
	                    if(polyverts[j][2]>Lz_-1){zplus=1;}
            		}
            	}

              double tempVec[3];
              double totVec[3]{0,0,0};
              double totArea = 0;
              double totalvert[3];
              double tempArea = 0;
              double cell0neighs[16]{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
              int cell0neighcount = 0;
              int counttemp = 0;

		    	for (long int i = 0; i < cell->polygons_.size(); i++) {  
	                
	                int numberOfVerticesOfFace = cell->polygons_[i]->vertices_.size();
	                double polyverts[numberOfVerticesOfFace][3];

                  tempVec[0]=0;
                  tempVec[1]=0;
                  tempVec[2]=0;
                  counttemp = 0;
                  tempArea = cell->polygons_[i]->area_;
                  totArea += cell->polygons_[i]->area_;
	                
	                //cout << "Number of Verticies of Polygon: " << polygons_[i]->vertices_.size() << endl;
	                for (int j = 0; j < cell->polygons_[i]->vertices_.size(); j++) {

	                    polyverts[j][0]=cell->polygons_[i]->vertices_[j]->position_[0];
	                    polyverts[j][1]=cell->polygons_[i]->vertices_[j]->position_[1];
	                    polyverts[j][2]=cell->polygons_[i]->vertices_[j]->position_[2];

	                    if(xminus==1 && xplus==1 && polyverts[j][0]<2){polyverts[j][0]+=Lx_;}
	                    if(yminus==1 && yplus==1 && polyverts[j][1]<2){polyverts[j][1]+=Ly_;}
	                    if(zminus==1 && zplus==1 && polyverts[j][2]<2){polyverts[j][2]+=Lz_;}

                      tempVec[0] += polyverts[j][0];
                      tempVec[1] += polyverts[j][1];
                      tempVec[2] += polyverts[j][2];
                      counttemp +=1;

	                    }             

                totVec[0] += tempVec[0]*tempArea/counttemp;
                totVec[1] += tempVec[1]*tempArea/counttemp;
                totVec[2] += tempVec[2]*tempArea/counttemp;

                double tempnorm[3];
                
                tempnorm[0]=0;
                tempnorm[1]=0;
                tempnorm[2]=0;
                int currpos;
                int nextpos;

                for (int j = 0; j < cell->polygons_[i]->vertices_.size(); j++) {
                	currpos = j;
                	nextpos = j+1;

                  if(j+1==cell->polygons_[i]->vertices_.size()){nextpos = 0;}

                	tempnorm[0] += (polyverts[currpos][1]-polyverts[nextpos][1])*(polyverts[currpos][2]+polyverts[nextpos][2]);
                	tempnorm[1] += (polyverts[currpos][2]-polyverts[nextpos][2])*(polyverts[currpos][0]+polyverts[nextpos][0]);
                	tempnorm[2] += (polyverts[currpos][0]-polyverts[nextpos][0])*(polyverts[currpos][1]+polyverts[nextpos][1]);

                }

                double tempmagvec;
                tempmagvec = pow(pow(tempnorm[0],2)+pow(tempnorm[1],2)+pow(tempnorm[2],2),0.5);
                cell->polygons_[i]->normvector_[0] = tempnorm[0]/tempmagvec;
                cell->polygons_[i]->normvector_[1] = tempnorm[1]/tempmagvec;
                cell->polygons_[i]->normvector_[2] = tempnorm[2]/tempmagvec;

              }
              

              totVec[0] = totVec[0]/totArea;
              totVec[1] = totVec[1]/totArea;
              totVec[2] = totVec[2]/totArea;

              for (long int i = 0; i < cell->polygons_.size(); i++) { 
              int curcellID = cell->id_;
              int curpolyID = cell->polygons_[i]->id_;

              if(cell->polygons_[i]->twincell_[0] ==-1){cell->polygons_[i]->twincell_[0] = cell->id_;
                  cell->polygons_[i]->twincell_[2] = curpolyID;}
              else if(cell->polygons_[i]->twincell_[2] == curpolyID){cell->polygons_[i]->twincell_[1] = cell->id_;
                  cell->polygons_[i]->twincell_[3] = curpolyID;}

              cell->polygons_[i]->cellcenter_[0] = totVec[0];
              cell->polygons_[i]->cellcenter_[1] = totVec[1];
              cell->polygons_[i]->cellcenter_[2] = totVec[2];

              cells_[cell->id_]->cellcenter_[0]  = totVec[0];
              cells_[cell->id_]->cellcenter_[1]  = totVec[1];
              cells_[cell->id_]->cellcenter_[2]  = totVec[2];      
         	 }
        }

        simulation_time_ += dt_;
    }


    return 0;
}

/*
 * Here, the moment of inertia considers a point mass at each vertex. Other formulations of
 * moment of inertia exist, but here we use the simplest which seems to work very well for our
 * needs of defining an ellipsoid for each cell. For more details, see Dobrovolskis (1996) - Inertia
 * of any polyhedron.
*/
EllipsoidByUnitPointMassPolyhedron Run::fitEllipsoidByUnitPointMassPolyhedron(int cellID) {
   
  // six terms of symmetric moment of intertia tensor
  double Ixx = 0;
  double Ixy = 0;
  double Ixz = 0;
  double Iyy = 0;
  double Iyz = 0;
  double Izz = 0;
  //int cellID = 1;


  int vertexCounter = 0;
  //for(DirectedFace *face : faces()) {
  for (long int i = 0; i < cells_[cellID]->polygons_.size(); i++) { 

    //DirectedEdgeOfCell *edge = face->firstEdge();
    //do {
    for (int j = 0; j < cells_[cellID]->polygons_[i]->vertices_.size(); j++) { 
      //edge = edge->nextAroundFace();
      //VertexOfCell *vertex = edge->vertex();

      double absx = cells_[cellID]->polygons_[i]->vertices_[j]->position_[0];
      double absy = cells_[cellID]->polygons_[i]->vertices_[j]->position_[1];
      double absz = cells_[cellID]->polygons_[i]->vertices_[j]->position_[2];
    
      double x = absx-cells_[cellID]->cellcenter_[0];
      double y = absy-cells_[cellID]->cellcenter_[1];
      double z = absz-cells_[cellID]->cellcenter_[2];

      if (x>Lx_/2) { x += -Lx_;}
      if (y>Ly_/2) { y += -Ly_;}
      if (z>Lz_/2) { z += -Lz_;}
      if (x<-Lx_/2) { x += Lx_;}
      if (y<-Ly_/2) { y += Ly_;}
      if (z<-Lz_/2) { z += Lz_;}

      Ixx += y*y + z*z;
      Iyy += x*x + z*z;
      Izz += x*x + y*y;
      Ixy += x*y;
      Ixz += x*z;
      Iyz += y*z;

      vertexCounter++;
      
      }
    //} while(edge!=face->firstEdge());
  } //end face loop
  
  Ixx /= vertexCounter;
  Iyy /= vertexCounter;
  Izz /= vertexCounter;
  Ixy /= vertexCounter;
  Ixz /= vertexCounter;
  Iyz /= vertexCounter;
  
  _unitPointMassMomentOfInertiaTensor = Matrix3x3(Vector3D(Ixx,-Ixy,-Ixz), Vector3D(-Ixy,Iyy,-Iyz), 
                                     Vector3D(-Ixz,-Iyz,Izz)); 
  
  return EllipsoidByUnitPointMassPolyhedron(_unitPointMassMomentOfInertiaTensor);
}


int     Run::updateVerticesVelocity() {
    for (auto vertex : vertices_) {
        for (int m = 0; m < 3; m++) {
            vertex->velocity_[m] = mu_ * (vertex->volumeForce_[m] + vertex->interfaceForce_[m]);
        }
    }
    //Account for spring verticies
    // for (long int i = 0; i < polygons_.size(); i++) {
    //     polygons_[i]->ResetTension();}
    for (long int i = 0; i < vertices_.size(); i++) {
        vertices_[i]->resetTensions();}
    for (long int i = 0; i < polygons_.size(); i++) {
        polygons_[i]->SpringGon(simulation_time_);}

    
    // remove drift velocity
    double averageVelocity[3] = {0., 0., 0.};
    for (auto vertex : vertices_) {
        for (int m = 0; m < 3; m++) {
            averageVelocity[m] = averageVelocity[m] + vertex->velocity_[m];
        }
    }
    for (int m = 0; m < 3; m++) {
        averageVelocity[m] = averageVelocity[m] / vertices_.size();
    }
    for (auto vertex : vertices_) {
        for (int m = 0; m < 3; m++) {
            vertex->velocity_[m] = vertex->velocity_[m] - averageVelocity[m];
        }
    }

    return 0;
}

int     Run::updateVerticesPosition() {
    std::default_random_engine generator(std::random_device{}());
    std::normal_distribution<double> ndist(0., 1.);
    std::uniform_real_distribution<> phinoise(-1, 1);
    double cR = sqrt(2.0*mu_*kB_*temperature_*dt_);
    double Dr = 1.0;
    double NoiseStdDev = sqrt(2*Dr*dt_);

    for (auto cell : cells_) {
        double tempvec[3];
        double tempvel[3];
        double tempnorm[2];

        tempvec[0] = ndist(generator);
        tempvec[1] = ndist(generator);
        tempvec[2] = ndist(generator);

        tempnorm[0] = sqrt(tempvec[0]*tempvec[0]+tempvec[1]*tempvec[1]+tempvec[2]*tempvec[2]);
        tempvec[0] = tempvec[0]/tempnorm[0];
        tempvec[1] = tempvec[1]/tempnorm[0];
        tempvec[2] = tempvec[2]/tempnorm[0];

        //if(cell->id_==1){cout << cell->cellDirectors_[0] << " " << cell->cellDirectors_[1] << endl;}

        tempvel[0] = cos(cell->cellDirectors_[1])*sin(cell->cellDirectors_[0]);
        tempvel[1] = sin(cell->cellDirectors_[1])*sin(cell->cellDirectors_[0]);
        tempvel[2] = cos(cell->cellDirectors_[0]);

        tempvel[0] += NoiseStdDev*(tempvec[0]-tempvel[0]);
        tempvel[1] += NoiseStdDev*(tempvec[1]-tempvel[1]);
        tempvel[2] += NoiseStdDev*(tempvec[2]-tempvel[2]);
        
        tempnorm[1] = sqrt(tempvel[0]*tempvel[0]+tempvel[1]*tempvel[1]+tempvel[2]*tempvel[2]);
        tempvel[0] = tempvel[0]/tempnorm[1];
        tempvel[1] = tempvel[1]/tempnorm[1];
        tempvel[2] = tempvel[2]/tempnorm[1];

        cell->cellDirectors_[0] = acos(tempvel[2]);
        cell->cellDirectors_[1] = atan2(tempvel[1],tempvel[0]);

        //if(cell->id_==1){cout << cell->cellDirectors_[0] << " " << cell->cellDirectors_[1] << endl;}

        //[Theta,Phi] Bias away from z
        //cell->cellDirectors_[0] += NoiseStdDev*(ndist(generator));
        //cell->cellDirectors_[1] += NoiseStdDev*(ndist(generator));

        //Bias in 
        //cell->cellDirectors_[0] += Dr/tan(cell->cellDirectors_[0]) + Randnormal(NoiseStdDev);
        //cell->cellDirectors_[1] += Randnormal(NoiseStdDev)/sin(cell->cellDirectors_[0]);
    }

     //for (long int i = 0; i < polygons_.size(); i++) {
     //    polygons_[i]->UpDateTension();}
     for (long int i = 0; i < cells_.size(); i++) {
        cells_[i]->updateAverageForce();}

    for (long int i = 0; i < vertices_.size(); i++) {
    	vertices_[i]->updateSP(temperature_,temperaturebot_,simulation_time_,placodeon_);
    	//vertices_[i]->updateVpos(simulation_time_);

        for (int m = 0; m < 3; m++) {
           //vertices_[i]->position_[m] = vertices_[i]->position_[m] + vertices_[i]->velocity_[m] * dt_ + cR*ndist(generator);
           vertices_[i]->position_[m] = vertices_[i]->position_[m] + vertices_[i]->velocity_[m] * dt_ + dt_ * vertices_[i]->motility_[m];
        }

        resetPosition(vertices_[i]->position_);
    }

    return 0;
}

double     Run::Randnormal(double stdev) {
    double rands = 0.0;
    rands = stdev*sqrt(-log(1.0-(rand()/((double)RAND_MAX+1)))*(2.0))*cos(2.0*M_PI*(rand()/((double)RAND_MAX+1)));

    return rands;
}

int     Run::updatePolygonVertices() {
    // update vertices in polygon
    for (long int i = 0; i < polygons_.size(); i++) {
        polygons_[i]->updateVertices();
    }

    return 0;
}

int     Run::updateCellVertices() {
    // update vertices in cell
    updatePolygonVertices();
    for (auto cell : cells_) {
        cell->vertices_.clear();
        for (auto polygon : cell->polygons_) {
            for (auto vertex : polygon->vertices_) {
                if (std::find(cell->vertices_.begin(), cell->vertices_.end(), vertex) == cell->vertices_.end()) {
                    // new vertex to be added
                    cell->vertices_.push_back(vertex);
                }
            }
        }
    }

    return 0;
}

int     Run::updateVertexEdges() {
    for (long int i = 0; i < vertices_.size(); i++) {
        vertices_[i]->edges_.clear();
    }
    for (long int i = 0; i < edges_.size(); i++) {
        edges_[i]->vertices_[0]->edges_.push_back(edges_[i]);
        edges_[i]->vertices_[1]->edges_.push_back(edges_[i]);
    }

    return 0;
}

int     Run::updateVertexCells() {
    for (auto vertex : vertices_) {
        vertex->cells_.clear();
    }
    std::vector<Cell *> tmp_cells = cells_;
    for (auto cell : tmp_cells) {
        for (auto polygon : cell->polygons_) {
            for (auto edge : polygon->edges_) {
                for (auto vertex : edge->vertices_) {
                    if (std::find(vertex->cells_.begin(), vertex->cells_.end(), cell) == vertex->cells_.end()) {
                        // new cell to be added
                        vertex->cells_.push_back(cell);
                    }
                }
            }
        }
    }

//    for (long int i = 0; i < vertices_.size(); i++) {
//        printf("%d\n", vertices_[i]->cells_.size());
//    }

    return 0;
}

int     Run::updatePolygonCells() {
    for (long int i = 0; i < polygons_.size(); i++) {
        polygons_[i]->cells_.clear();
    }
    for (long int i = 0; i < cells_.size(); i++) {
        for (int j = 0; j < cells_[i]->polygons_.size(); j++) {
            cells_[i]->polygons_[j]->cells_.push_back(cells_[i]);
        }
    }

    return 0;
}

int     Run::updateCellShapeIndex() {
    for (auto cell : cells_) {
        double area = 0.;
        for (auto polygon : cell->polygons_) {
            area += polygon->area_;
        }
        cell->shapeIndex_ = area * pow(cell->volume_, (-1.0)*2.0/3.0);
    }

    return 0;
}

int     Run::updateGeoinfo() {
    // update edge midpoint and length
    for (long int i = 0; i < edges_.size(); i++) {
        edges_[i]->update();
//        printf("%6f\n", run->edges_[i]->length_);
    }
    // update polygon center position
    for (long int i = 0; i < polygons_.size(); i++) {
        polygons_[i]->updateCenter();
    }

    return 0;
}

int     Run::deleteVertex(Vertex * vertex) {
    auto it = find(vertices_.begin(), vertices_.end(), vertex);
    if (it != vertices_.end()) {
//        int index = it - vertices_.begin();
        vertices_.erase(it);
    } else {
        printf("vertex %ld not found in vertices_\n", vertex->id_);
        exit(1);
    }
    delete vertex;

    return 0;
}

int     Run::deleteEdge(Edge * edge) {
    auto it = find(edges_.begin(), edges_.end(), edge);
    if (it != edges_.end()) {
        edges_[it-edges_.begin()]->markToDelete_ = true;
//        edges_.erase(it);
    } else {
        printf("edge %ld not found in edges_\n", edge->id_);
        exit(1);
    }
//    delete edge;

    return 0;
}

int     Run::deletePolygon(Polygon * polygon) {
    auto it = find(polygons_.begin(), polygons_.end(), polygon);
    if (it != polygons_.end()) {
//        int index = it - vertices_.begin();
        polygons_.erase(it);
    } else {
        printf("polygon %ld not found in polygons_\n", polygon->id_);
        exit(1);
    }
    delete polygon;

    return 0;
}

int     Run::resetPosition(double * r) {
    if (fabs(r[0]) > 1e6) {
        printf("%e\n",r[0]);
    }
    if (fabs(r[1]) > 1e6) {
        printf("%e\n",r[1]);
    }
    if (fabs(r[2]) > 1e6) {
        printf("%e\n",r[2]);
    }
    r[0] = r[0] - Lx_ * floor(r[0] / Lx_);
    r[1] = r[1] - Ly_ * floor(r[1] / Ly_);
    r[2] = r[2] - Lz_ * floor(r[2] / Lz_);

    return 0;
}

Edge *  Run::addEdge(Vertex * v0, Vertex * v1) {
    Edge * edge = new Edge(this, count_edges_);
    count_edges_ += 1;
    edges_.push_back(edge);
    if (v0->id_ < v1->id_) {
        edge->vertices_.push_back(v0);
        edge->vertices_.push_back(v1);
    } else {
        edge->vertices_.push_back(v1);
        edge->vertices_.push_back(v0);
    }
    edge->update();
    edge->candidate_ = false;

    return edge;
}

int Run::dumpConfigurationVtk() {
    //////////////////////////////////////////////////////////////////////////////////////
    stringstream filename;
    filename << setw(7) << setfill('0') << (long int)(floor(simulation_time_+0.01*dt_)) << ".sample.vtk";
    ofstream out(filename.str().c_str());
    if (!out.is_open()) {
        cout << "Error opening output file " << filename.str().c_str() << endl;
        exit(1);
    }
    out << "# vtk DataFile Version 2.0" << endl;
    out << "polydata" << endl;
    out << "ASCII" << endl;
    out << "DATASET POLYDATA" << endl;
    //out << "DATASET UNSTRUCTURED_GRID" << endl;

    out << "POINTS " << vertices_.size() << " double" << endl;
    for (long int i = 0; i < vertices_.size(); i++) {
        // reset vertex id for dumping polygons
        vertices_[i]->dumpID_ = i;
        out << right << setw(12) << scientific << setprecision(5) << vertices_[i]->position_[0];
        out << " " << right << setw(12) << scientific << setprecision(5) << vertices_[i]->position_[1];
        out << " " << right << setw(12) << scientific << setprecision(5) << vertices_[i]->position_[2];
        out << endl;
    }
    out << endl;


//    long int Nedges = 0;
//    for (long int i = 0; i < run->edges_.size(); i++) {
//        if (!run->edges_[i]->crossBoundary()) {
//            Nedges++;
//        }
//    }
//    out << "LINES " << Nedges << " " << 3*Nedges << endl;
//    for (long int i = 0; i < run->edges_.size(); i++) {
//        if (!run->edges_[i]->crossBoundary()) {
//            out << left << setw(6) << 2;
//            for (int j = 0; j < run->edges_[i]->vertices_.size(); j++) {
//                out << " " << left << setw(6) << run->edges_[i]->vertices_[j]->id_;
//            }
//            out << endl;
//        }
//    }
//    out << endl;

    updatePolygonVertices();
    long int Npolygons = 0;
    long int NpolygonVertices = 0;
    for (long int i = 0; i < polygons_.size(); i++) {
        if (!polygons_[i]->crossBoundary()) {
            Npolygons++;
            NpolygonVertices += polygons_[i]->vertices_.size();
        }
    }
    out << "POLYGONS " << Npolygons << " " << Npolygons + NpolygonVertices << endl;
    for (long int i = 0; i < polygons_.size(); i++) {
        if (!polygons_[i]->crossBoundary()) {
            out << left << setw(6) << polygons_[i]->vertices_.size();
            for (int j = 0; j < polygons_[i]->vertices_.size(); j++) {
                out << " " << left << setw(6) << polygons_[i]->vertices_[j]->dumpID_;
            }
            out << endl;
        }
    }
    out << endl;

    // out << "CELL_TYPES " << Npolygons << endl;
    // for (long int i = 0; i < polygons_.size(); i++) {
    //     out << left << setw(6) << polygons_[i]->type_;
    //     out << endl;
    // }
    // out << endl;

//    updatePolygonVolumeRatio();
//    out << "CELL_DATA " << Npolygons << endl;
//    out << "SCALARS volumeRatio double 1" << endl;
//    out << "LOOKUP_TABLE default" << endl;
//    for (long int i = 0; i < polygons_.size(); i++) {
//        if (!polygons_[i]->crossBoundary()) {
//            out << left << setw(6) << polygons_[i]->dumpVolumeRatio_ << endl;
//        }
//    }
//    out << endl;

    // updateCellVertices();
    // long int NCells = 0;
    // long int NCellVerts = 0;
    // for (long int i = 0; i < cells_.size(); i++) {
    //     NCells++;
    //     NCellVerts += cells_[i]->vertices_.size();

    // }
    // out << "CELL_DATA  " << NCells << " " << NCells + NCellVerts << endl;
    // for (long int i = 0; i < cells_.size(); i++) {
    //     out << left << setw(6) << cells_[i]->vertices_.size();
    //     for (int j = 0; j < cells_[i]->vertices_.size(); j++) {
    //         out << " " << left << setw(6) << cells_[i]->vertices_[j]->dumpID_;
    //     }
    //     out << endl;
    // }
    // out << endl;

    // out << "CELL_TYPES " << NCells << endl;
    // for (long int i = 0; i < cells_.size(); i++) {
    //     out << left << setw(6) << cells_[i]->type_;
    //     out << endl;
    // }
    // out << endl;



    out.close();

    return 0;
}

int     Run::dumpCellCenter() {
    updateCellVertices();
    stringstream filename;
    filename << "cellCenter.txt";
    ofstream out(filename.str().c_str(), std::ios_base::app);
    if (!out.is_open()) {
        cout << "Error opening output file " << filename.str().c_str() << endl;
        exit(1);
    }
    out << "time ";
    out << left << setw(12) << simulation_time_;
    out << endl;

    for (auto cell : cells_) {
        double center[3] = {0., 0., 0.};
        double reference[3];
        for (int m = 0; m < 3; m++) {
            reference[m] = cell->vertices_[0]->position_[m];
        }
        for (auto vertex : cell->vertices_) {
            double dx[3];
            for (int m = 0; m < 3; m++) {
                dx[m] = (vertex->position_[m] - reference[m]);
            }
            while (dx[0] > Lx_/2.0) {
                dx[0] = dx[0] - Lx_;
            }
            while (dx[0] < (-1.0)*Lx_/2.0) {
                dx[0] = dx[0] + Lx_;
            }
            while (dx[1] > Ly_/2.0) {
                dx[1] = dx[1] - Ly_;
            }
            while (dx[1] < (-1.0)*Ly_/2.0) {
                dx[1] = dx[1] + Ly_;
            }
            while (dx[2] > Lz_/2.0) {
                dx[2] = dx[2] - Lz_;
            }
            while (dx[2] < (-1.0)*Lz_/2.0) {
                dx[2] = dx[2] + Lz_;
            }
            for (int m = 0; m < 3; m++) {
                center[m] = center[m] + dx[m];
            }
        }
        for (int m = 0; m < 3; m++) {
            center[m] = center[m]/cell->vertices_.size() + reference[m];
        }
        resetPosition(center);
        out << left << setw(6) << cell->id_;
        out << " " << right << setw(12) << scientific << setprecision(5) << center[0];
        out << " " << right << setw(12) << scientific << setprecision(5) << center[1];
        out << " " << right << setw(12) << scientific << setprecision(5) << center[2];
        out << endl;
    }
    out << endl;

    out.close();

    return 0;
}

int     Run::dumpCellShapeIndex() {
    updateCellShapeIndex();
    stringstream filename;
    filename << "cellShapeIndex.txt";
    ofstream out(filename.str().c_str(), std::ios_base::app);
    if (!out.is_open()) {
        cout << "Error opening output file " << filename.str().c_str() << endl;
        exit(1);
    }
    out << "time ";
    out << left << setw(12) << simulation_time_;
    double averageShapeIndex = 0.;
    for (auto cell : cells_) {
        averageShapeIndex += cell->shapeIndex_;
    }
    averageShapeIndex /= cells_.size();
    out << " " << left << setw(12) << averageShapeIndex;
    out << endl;

    for (auto cell : cells_) {
        out << left << setw(6) << cell->id_;
        out << " " << cell->shapeIndex_;
        out << endl;
    }
    out << endl;

    out.close();

    return 0;
}


int     Run::dumpDemix() {
    stringstream filename;
    filename << "demixing.txt";
//    ofstream out(filename.str().c_str(), std::ios::binary | std::ios_base::app);
    ofstream out(filename.str().c_str(), std::ios_base::app);
    if (!out.is_open()) {
        cout << "Error opening output file " << filename.str().c_str() << endl;
        exit(1);
    }

    out << "time ";
    out << left << setw(12) << simulation_time_;
    double averageShapeIndex = 0.;

    double demixing = 0.;
    double countcells = 0.;
    for (auto cell : cells_) 
        {
        countcells += 1;
        double countpoly = 0.;
        double counthom = 0.;
        for (auto polygon : cell->polygons_) 
            {
            countpoly += 1;
            if(polygon->type_!=1)
                {
                    counthom+=1;
                    //cout << polygon->type_ << endl;
                }   
            }
        demixing += 2.0*((counthom/countpoly)-0.5);
        //cout << ((counthom/countpoly)-0.5) << endl;
        }

    out << "demix ";
    out << left << setw(12) << demixing/(countcells);
    out << endl;

    out << endl;
    out.close();

    return 0;
}


int     Run::dumpTopo() {
    stringstream filename;
    filename << "topo.txt";
//    ofstream out(filename.str().c_str(), std::ios::binary | std::ios_base::app);
    ofstream out(filename.str().c_str(), std::ios_base::app);
    if (!out.is_open()) {
        cout << "Error opening output file " << filename.str().c_str() << endl;
        exit(1);
    }
    out << "time ";
    out << left << setw(12) << simulation_time_;
    out << endl;

    out << "vertices ";
    out << left << setw(12) << vertices_.size();
    out << endl;
    for (auto vertex : vertices_) {
        out << left << setw(6) << vertex->id_;
        out << " " << right << setw(12) << scientific << setprecision(5) << vertex->position_[0];
        out << " " << right << setw(12) << scientific << setprecision(5) << vertex->position_[1];
        out << " " << right << setw(12) << scientific << setprecision(5) << vertex->position_[2];
        out << endl;
    }

    out << "edges ";
    out << left << setw(12) << edges_.size();
    out << endl;
    for (auto edge : edges_) {
        out << left << setw(6) << edge->id_;
        for (auto vertex : edge->vertices_) {
            out << " " << right << setw(12) << scientific << setprecision(5) << vertex->id_;
        }
        out << endl;
    }

    out << "polygons ";
    out << left << setw(12) << polygons_.size();
    out << endl;
    for (auto polygon : polygons_) {
        out << left << setw(6) << polygon->id_;
        for (auto edge : polygon->edges_) {
            out << " " << right << setw(12) << scientific << setprecision(5) << edge->id_;
        }
        out << endl;
    }

    out << "cells ";
    out << left << setw(12) << cells_.size();
    out << endl;
    for (auto cell : cells_) {
        out << left << setw(6) << cell->id_;
        for (auto polygon : cell->polygons_) {
            out << " " << right << setw(12) << scientific << setprecision(5) << polygon->id_;
        }
        out << endl;
    }

    out << endl;
    out.close();

    return 0;
}

int     Run::dumpReconnection() {
    stringstream filename;
    filename << "reconnections.txt";
    ofstream out(filename.str().c_str(), std::ios_base::app);
    if (!out.is_open()) {
        cout << "Error opening output file " << filename.str().c_str() << endl;
        exit(1);
    }

    out << verboseReconnection_.str();
    verboseReconnection_.str("");

    out.close();

    return 0;
}
