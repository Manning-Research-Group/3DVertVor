/* 
 * File:   MovieCreator.h
 * Author: arbeit
 *
 * Created on October 29, 2015, 5:25 PM
 */

#ifndef MOVIECREATOR_H
#define	MOVIECREATOR_H

#include <string>

class MovieCreator {
public:
  MovieCreator();
  ~MovieCreator();
  bool start();
  std::string nextFrameFileName();
  bool createMovie(const std::string &fileName, int rate=30);
  bool cleanUp();
  
private:
  const static char BaseFolder[];
  const static char FolderMask[];
  const static char FrameMask[];
  const static char LogFileName[];
  const static int MaxPathNameLength;

  bool _valid;
  std::string _folderName;
  std::string _fileNameMask;
  int _frameCounter;
};

#endif	/* MOVIECREATOR_H */

