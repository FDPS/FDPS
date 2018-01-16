#ifndef H_PDBMANAGER
#define H_PDBMANAGER

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <cstdio>

class PDBData{
public:
  std::string  record;
  //12345678901234567890
  //
  unsigned int serial;   // atom serial number
  std::string  name;     // atom name
  std::string  alt_loc;  // alternate location indicator
  std::string  resName;  // residue name
  std::string  chainID;  // chain indicator
  unsigned int resSeq;   // residue sequential number
  std::string  iCode;    // code for insersion of residues
  double       coord[3]; // orthogonal coordinate in angstroms
  double       occupancy;// occupancy
  double       tempFac;  // temperature factor
  std::string  element;  // element symbol
  std::string  charge;   // charge
  PDBData(){};
  PDBData(const PDBData& _pdb){
    serial    = _pdb.serial;   
    name      = _pdb.name;     
    alt_loc   = _pdb.alt_loc;  
    resName   = _pdb.resName;  
    chainID   = _pdb.chainID;  
    resSeq    = _pdb.resSeq;   
    iCode     = _pdb.iCode;    
    coord[0]  = _pdb.coord[0];
    coord[1]  = _pdb.coord[1];
    coord[2]  = _pdb.coord[2]; 
    occupancy = _pdb.occupancy;
    tempFac   = _pdb.tempFac;  
    element   = _pdb.element;  
    charge    = _pdb.charge;   
  }
  
  friend std::istream& operator>>(std::istream &is,PDBData &d){
    std::string line;
    std::getline(is,line);
    if(line == "END"){
      d.record = line;
      return is;
    }

    d.record = line.substr(0,6);
    if(d.record == "REMARK") return is;
    if(d.record=="HETATM" || d.record == "ATOM"){
      d.serial    = std::stoi(line.substr( 6,5));
      d.name      =           line.substr(12,4);
      d.alt_loc   =           line.substr(16,1);
      d.resName   =           line.substr(17,3);
      d.chainID   =           line.substr(21,1);
      d.resSeq    = std::stoi(line.substr(22,4));
      d.iCode     =           line.substr(26,1);
      d.coord[0]  = std::stod(line.substr(30,8));
      d.coord[1]  = std::stod(line.substr(38,8));
      d.coord[2]  = std::stod(line.substr(46,8));
      d.occupancy = std::stod(line.substr(54,6));
      d.tempFac   = std::stod(line.substr(60,6));
      d.element   =           line.substr(76,2);
      d.charge    =           line.substr(78,2);
      //std::cout << d << std::endl;
      return is;
    }
    assert(false);
  }
  friend std::ostream& operator<<(std::ostream &os,PDBData &d){
    os << std::left << std::setw(6); os << d.record;
    os << std::right;
    if(d.record=="HETATM"){
      os << std::setw(5); os << d.serial;
      os << " ";
      os << std::setw(4); os << d.name;
      os << std::setw(1); os << d.alt_loc;
      os << std::setw(3); os << d.resName;
      os << " ";
      os << std::setw(1); os << d.chainID;
      os << std::setw(4); os << d.resSeq;
      os << std::setw(1); os << d.iCode;
      os << "   ";
      os << std::fixed << std::setprecision(3);
      os << std::setw(8); os << d.coord[0];
      os << std::setw(8); os << d.coord[1];
      os << std::setw(8); os << d.coord[2];
      os << std::fixed << std::setprecision(2);
      os << std::setw(6); os << d.occupancy;
      os << std::setw(6); os << d.tempFac;
      os << "          ";
      os << std::setw(2); os << d.element;
      os << std::setw(2); os << d.charge;
    }
    return os;
  }
};

template <class Data>
class FileManager{
private:
public:
  Data ReadLine(const std::string line){
    Data ret;
    if(line == "END" || line=="TER"){
      ret.record = line;
      return ret;
    }
    ret.record = line.substr(0,6);
    if(ret.record == "REMARK") return ret;
    if(ret.record=="HETATM" || ret.record == "ATOM"){
      ret.serial    = std::stoi(line.substr( 6,5));
      ret.name      =           line.substr(12,4);
      ret.alt_loc   =           line.substr(16,1);
      ret.resName   =           line.substr(17,3);
      ret.chainID   =           line.substr(21,1);
      ret.resSeq    = std::stoi(line.substr(22,4));
      ret.iCode     =           line.substr(26,1);
      ret.coord[0]  = std::stod(line.substr(30,8));
      ret.coord[1]  = std::stod(line.substr(38,8));
      ret.coord[2]  = std::stod(line.substr(46,8));
      ret.occupancy = std::stod(line.substr(54,6));
      ret.tempFac   = std::stod(line.substr(60,6));
      ret.element   =           line.substr(76,2);
      ret.charge    =           line.substr(78,2);
      //std::cout << d << std::endl;
      return ret;
    }
    std::cerr << "error: unsupported PDB record type!" << std::endl;
    assert(false);
  }
  void WriteLine(std::istream &is){
    Data ret; is << ret;
    return ret;
  }
};

#endif
