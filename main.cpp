#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <unordered_set>
#include <filesystem>
#include <fstream>
#include "MurmurHash3.h"

using namespace std;
using namespace std::filesystem;

const string g_path = "/home/ralsei/Escritorio/Topicos en Volumen de Datos/Trabajo1/Genomes/";

class HyperLogLog{
    private: 
        int b;
        int m;
        vector<int>registro;
        double alpha()const{
            switch(this->m){
                case 16: return 0.673;
                case 32: return 0.697;
                case 64: return 0.709;
                default: return 0.7213/(1+1.079/this->m);
            }
        }

    public:
        HyperLogLog(int b){
            this->b = b;
            this->m = pow(2,b);
            registro.resize(this->m,0);
        }
        ~HyperLogLog(){
            registro.clear();
        }
        void ingresar(string s){
            uint32_t hash;
            MurmurHash3_x86_32(s.c_str(), s.size(), 0, &hash);
            int indice = hash >> (32 - b);
            int rango = __builtin_clz((hash << b) | (1 << (b - 1))) + 1;
            registro[indice] = max(registro[indice], rango);

            // Verificar si el string se ingresó correctamente
            //cout << "String: " << s << " | Hash: " << hash << " | Indice: " << indice << " | Rango: " << rango << " | Registro[indice]: " << registro[indice] << endl;
        }
        double estimar(){
            double alfaMM = alpha()*this->m*this->m;
            double suma = 0.0;
            for(int i = 0; i < this->m; i++){
                suma += 1.0/(1 << registro[i]);
            }
            double estimacion = alfaMM/suma;
            //cout << "Estimacion pre-correccion: " << estimacion << endl;
            if(estimacion <= 5.0/2.0*this->m){
                double ceros = 0;
                for(int i = 0; i < this->m; i++){
                    if(registro[i] == 0){
                        ceros++;
                    }
                }
                if(ceros != 0){
                    //cout << "Correccion: " << this->m << " * log(" << this->m << "/" << ceros << ")" << endl;
                    //cout << this->m * log(this->m/ceros) << endl;
                    return this->m*log(this->m/ceros);
                }
            }else if(estimacion > 1.0/30.0*pow(2,32)){
                //cout << "Correccion: " << -pow(2,32) << " * log(1 - " << estimacion << "/" << pow(2,32) << ")" << endl;
                return -pow(2,32)*log(1-estimacion/pow(2,32));
            }
            return estimacion;
        }
        vector<int> getRegistro(){
            return registro;
        }
        void combinar(HyperLogLog hll1, HyperLogLog hll2){
            vector<int> registro1 = hll1.getRegistro();
            vector<int> registro2 = hll2.getRegistro();
            //cout << "Combinando" << endl;
            for(int i = 0; i < this->m; i++){
                this->registro[i] = max(registro1[i], registro2[i]);
            }
            //cout << "Combinacion exitosa" << endl;
        }
};

vector<string> K_merizadorInador(const string &s, int k){
    vector<string> k_mers;
    for(int i = 0; i < s.size() - k + 1; i++){
        k_mers.push_back(s.substr(i,k));
    }
    return k_mers;
}

vector<string> MinimizadorInador(const string &s, int k, int w){
    vector<string> minimizadores;
    for(int i=0; i < s.size()-w+1; i++){
        string subcadena = s.substr(i,w);
        //cout << "Subcadena: " << subcadena << endl;
        vector<string> k_mers = K_merizadorInador(subcadena,k);
        /*cout << "K_mers candidatos: " << k_mers.size() << endl;
        for(auto &k_mer: k_mers){
            cout << "K_mer: " << k_mer << endl;
        }*/
        sort(k_mers.begin(), k_mers.end());
        minimizadores.push_back(k_mers[0]);
    }
    return minimizadores;
}

double jaccard(HyperLogLog hll1, HyperLogLog hll2, int b){
    double estimacion1 = hll1.estimar();
    //cout << "Estimacion1: " << estimacion1 << endl;
    double estimacion2 = hll2.estimar();
    //cout << "Estimacion2: " << estimacion2 << endl;

    HyperLogLog hll_union(b);
    hll_union.combinar(hll1,hll2);
    double estimacion_union = hll_union.estimar();
    //cout << "Estimacion_union: " << estimacion_union << endl;
    return (estimacion1 + estimacion2 - estimacion_union)/estimacion_union;
}

int main(int argc, char *argv[]){
    int k = 20;
    int w = 50;
    int b = 16;
    /*string test = "ACGTGACCG";
    vector<string> minimizadores = MinimizadorInador(test, 3, 6);
    for(auto &minimizador: minimizadores){
        cout << minimizador << endl;
    }*/
    vector<string> genomas;
    for(auto &p: directory_iterator(g_path)){
        ifstream file(p.path());
        string s;
        string genoma = "";
        while(getline(file,s)){
            if(s[0] == '>'){
                continue;
            }
            genoma += s;
        }
        genomas.push_back(genoma);
        
    }
    vector<HyperLogLog> hlls;
    vector<unordered_set<string>> exactos;
    for(int i = 0; i < genomas.size(); i++){
        HyperLogLog hll(b);
        vector<string> kmers1 = K_merizadorInador(genomas[i], k);
        vector<string> minimizadores1 = MinimizadorInador(genomas[i], k, w);
        cout << "Kmers: " << kmers1.size() << " Minimizadores: " << minimizadores1.size() << endl;
        for(auto &kmer: kmers1){
            hll.ingresar(kmer);
        }
        HyperLogLog hll2(b);
        for(auto &minimizador: minimizadores1){
            hll2.ingresar(minimizador);
        }
        hlls.push_back(hll);
        hlls.push_back(hll2);
    }

    for(int i = 0; i < hlls.size(); i+=2){
        if(i+2 >= hlls.size()){
            break;
        }
        cout << "Jaccard kmers: " << jaccard(hlls[i],hlls[i+2],b) << endl;
        cout << "Jaccard minimizadores: " << jaccard(hlls[i+1],hlls[i+3],b) << endl;
    }
    hlls.clear();
    for(int i = 0; i < genomas.size(); i++){
        unordered_set<string> exacto;
        vector<string> kmers = K_merizadorInador(genomas[i], k);
        vector<string> minimizadores = MinimizadorInador(genomas[i], k, w);
        for(auto &kmer: kmers){
            exacto.insert(kmer);
        }
        exactos.push_back(exacto);
        exacto.clear();
        for(auto &minimizador: minimizadores){
            exacto.insert(minimizador);
        }
        exactos.push_back(exacto);
    }
    for(int i = 0; i < exactos.size(); i+=2){
        if(i+2 >= exactos.size()){
            break;
        }
        unordered_set<string> interseccion;
        for(const auto& elem : exactos[i]) {
            if (exactos[i+2].find(elem) != exactos[i+2].end()) {
            interseccion.insert(elem);
            }
        }
        cout << "Jaccard exactos kmers: " << (double)interseccion.size()/(exactos[i].size() + exactos[i+2].size() - interseccion.size()) << endl;
        interseccion.clear();
        for(const auto& elem : exactos[i+1]) {
            if (exactos[i+3].find(elem) != exactos[i+3].end()) {
            interseccion.insert(elem);
            }
        }
        cout << "Jaccard exactos minimizadores: " << (double)interseccion.size()/(exactos[i+1].size() + exactos[i+3].size() - interseccion.size()) << endl;
    }
    
    return 0;
}