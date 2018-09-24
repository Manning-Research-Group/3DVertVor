/* 
 * File:   MovieCreator.cpp
 * Author: arbeit
 * 
 * Created on October 29, 2015, 5:25 PM
 */
#include <stdlib.h>
#include <sstream>

#include "fileIo.h"
#include <iostream>
#include "MovieCreator.h"

const char MovieCreator::BaseFolder[] = "temp";
const char MovieCreator::FolderMask[] = "temp/XXXXXX";
const char MovieCreator::FrameMask[] = "frame%04d.png";
const char MovieCreator::LogFileName[] = "ffmpeg.log";
const int MovieCreator::MaxPathNameLength = 200;

MovieCreator::MovieCreator() : _valid(false) {}
MovieCreator::~MovieCreator() { 
  cleanUp();
}

bool MovieCreator::start() {
  cleanUp();
  
  // create base folder
  if(!recursivelyCreateFolder(BaseFolder)) {
    return false;
  }
  
  // generate temporary path
  char pathName[MaxPathNameLength];
  strcpy(pathName, FolderMask);
  if(!mkdtemp(pathName)) {
    return false;
  }

  _folderName = pathName;
  _frameCounter = 0;
  
  char temp[MaxPathNameLength];
  sprintf(temp, "%s/%s", _folderName.c_str(), FrameMask);
  std::stringstream ss1;
  ss1 << _folderName << "/" << FrameMask;
  _fileNameMask = ss1.str();

  _valid = true;
  return true;
}

std::string MovieCreator::nextFrameFileName() {
  if(_valid) {
    char temp[MaxPathNameLength];
    sprintf(temp, _fileNameMask.c_str(), _frameCounter);
    ++_frameCounter;
    return temp;
  }
  return std::string();
}

bool MovieCreator::createMovie(const std::string &fileName, int rate) {
  if(_valid && (_frameCounter>0)) {
    std::stringstream ss;
    ss << "ffmpeg -r " << rate << " -i " << _fileNameMask << " -preset slow -crf 18 -f mp4 -vcodec h264 -pix_fmt yuv420p -y " << fileName << " 2>" << _folderName << "/" << LogFileName;
//    ss << "ffmpeg -r " << rate << " -i " << _fileNameMask << " -q:v 1 -y " << fileName << " 2>" << _folderName << "/" << LogFileName;
    std::cout << ss.str() << std::endl;
    if(0==system(ss.str().c_str())) {
//      std::cout << "Done." << std::endl;
      return true;
    } else {
      std::cout << "Error creating movie!" << std::endl;
      return false;
    }  
  }
  return false;
}


bool MovieCreator::cleanUp() {
  if(_valid) {
    recursivelyDeleteFolder(_folderName);
    _valid = false;
  }
  return true;
}


