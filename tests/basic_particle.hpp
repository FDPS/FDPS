#pragma once

class BasicParticle32 {
public:
    PS::S32    id;
    PS::F32vec pos;

    PS::F32vec getPos() const {
        return this->pos;
    }

    void generateSphere(PS::S64 i,
                        PS::F32 radius,
                        PS::F32vec center = 0.d) {
        this->id = i;
        do {
            for(PS::S32 k = 0; k < PS::DIMENSION; k++)
                this->pos[k] = radius * (2.0d * PS::MT::genrand_real2() - 1.0d);
        }while(this->pos * this->pos >= radius * radius);
        this->pos += center;
        this->pos += center;
    }

    void writeAscii(FILE *fp) {
        fprintf(fp, "%6d %+e %+e %+e\n", this->id, this->pos[0], this->pos[1], this->pos[2]);
    }

};

class BasicParticle64 {
public:
    PS::S64    id;
    PS::F64vec pos;

    PS::F64vec getPos() const {
        return this->pos;
    }

    void generateSphere(PS::S64 i,
                        PS::F64 radius,
                        PS::F64vec center = 0.d) {
        this->id = i;
        do {
            for(PS::S32 k = 0; k < PS::DIMENSION; k++)
                this->pos[k] = radius * (2.0d * PS::MT::genrand_real2() - 1.0d);
        }while(this->pos * this->pos >= radius * radius);
        this->pos += center;
    }

    void writeAscii(FILE *fp) {
        fprintf(fp, "%6d %+e %+e %+e\n", this->id, this->pos[0], this->pos[1], this->pos[2]);
    }

};
