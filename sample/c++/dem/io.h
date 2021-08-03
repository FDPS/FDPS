#pragma pnce

template <class ThisPtcl> class FileIO{
	PS::F64 time;
	PS::S64 step;
	public:
	FileIO(): time(0.0), step(0){
	}
	void OutputFileWithTimeInterval(PS::ParticleSystem<ThisPtcl>& ptcl, const system_t& sysinfo){
		const int NumberOfSnapshot = sysinfo.end_time * 60;//60 frames per sec
		if(sysinfo.time >= time){
			FileHeader header;
			header.time = sysinfo.time;
			header.Nbody = ptcl.getNumberOfParticleLocal();
			char filename[256];
			sprintf(filename, "result/%05d", step);
			ptcl.writeParticleAscii(filename, "%s_%05d_%05d.dat", header);
			if(PS::Comm::getRank() == 0){
				std::cerr << "output " << filename << "." << std::endl;
			}
			//std::cout << time << " " << ptcl[0].pos << " " << ptcl[0].vel << " " << ptcl[0].avel << std::endl;
			time += sysinfo.end_time / NumberOfSnapshot;
			++ step;
		}
	};
	//Memento
	bool Create(const PS::ParticleSystem<ThisPtcl>& ptcl, const system_t& sysinfo){
		std::ofstream fout;
		char filename[256];
		sprintf(filename, "result/%05d_%05d_%05d.bin", 0, PS::Comm::getNumberOfProc(), PS::Comm::getRank());
		fout.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
		if(!fout){
			std::cerr << "cannot open restart file." << std::endl;
			exit(1);
		}
		fout.write(reinterpret_cast<const char * const>(this), sizeof(FileIO));
		fout.write(reinterpret_cast<const char * const>(&sysinfo), sizeof(system_t));
		for(std::size_t i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
			const ThisPtcl& ith = ptcl[i];
			fout.write((char*)&ith, sizeof(ThisPtcl));
		}
		fout.close();
		std::cerr << "created Memento file" << std::endl;
		return true;
	}
	bool Restore(PS::ParticleSystem<ThisPtcl>& ptcl, system_t& sysinfo){
		char filename[256];
		sprintf(filename, "result/%05d_%05d_%05d.bin", 0, PS::Comm::getNumberOfProc(), PS::Comm::getRank());
		std::ifstream fin(filename, std::ios::in | std::ios::binary);
		if(!fin){
			std::cerr << "cannot open restart file." << std::endl;
			exit(1);
		}
		std::vector<ThisPtcl> ptcl_loc;
		fin.read((char*)this, sizeof(FileIO));
		fin.read((char*)&sysinfo, sizeof(system_t));
		while(1){
			ThisPtcl ith;
			fin.read((char*)&ith, sizeof(ThisPtcl));
			if(fin.eof() == true) break;
			ptcl_loc.push_back(ith);
		}
		fin.close();
		ptcl.setNumberOfParticleLocal(ptcl_loc.size());
		for(std::size_t i = 0 ; i < ptcl_loc.size() ; ++ i){
			ptcl[i] = ptcl_loc[i];
		}
		std::cerr << "IP" << std::endl;
		return true;
	}
};

