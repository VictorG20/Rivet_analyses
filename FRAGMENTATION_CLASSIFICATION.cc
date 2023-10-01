// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <filesystem>  // Used to create directory.
#include <iomanip>  // std::setprecision



#include "Rivet/Event.hh"
#include "Rivet/Particle.hh"

namespace Rivet {

/// @brief Underlying event in Z events
class FRAGMENTATION_CLASSIFICATION : public Analysis {
public:
    
    inline bool exists (const std::string& filename) {
      struct stat buffer;
      return (stat (filename.c_str(), &buffer) == 0);
    }
    
    template <typename T>
    void write_vector(ofstream& outfile, string header, const std::vector<T> &vect) {
        
        outfile << std::setw(my_width1) << std::left << header << ", ";
        
        for (int i = 0; i < vect.size(); i++){
            
            outfile << std::setprecision(my_width2 - 3);
            
            if (i == vect.size() - 1) {
                outfile << std::setw(my_width2) << std::right << vect[i] << "\n";
                continue;
            }
            
            outfile << std::setw(my_width2) << std::right << vect[i] << ", ";
        }
        
        return;
    }
    
    // -----------------------------------------------------------------------------------------------------------------
    // Function to classify the particles incoming to the Fragmentation vertex according to
    // whether or not they are children to the outgoing particles of the signal process.
    // -----------------------------------------------------------------------------------------------------------------
    void classify_particles(const Event& event, ofstream& outfile, double Zphi, double Zeta, double Zpt) {
        int sign_status = 1;  // Signal process vertex has status = 1.
        int hcol_status = 2;  // Hard collision vertices have status = 2.
        int frag_status = 5;  // Fragmentation vertex has status = 5.
        int hcol_counter = 0;  // Hard collisions counter.
        
        std::vector<int> taken_ids, signals, pid_values;
        std::vector<float> deta_values, dphi_values, pT_values;
        
        double pTsumTow(0.0), pTsumTrans(0.0), pTsumAway(0.0);
        double pTsumTransmin(0.0), pTsumTransmax(0.0), pTsumLeft(0.0), pTsumRight(0.0);
        
        for(auto & vertex : (*event.genEvent()).vertices()){
            int v_stat = vertex -> status();
            
            // Signal process vertex always comes before the Fragmentation vertex.
            if(v_stat == sign_status){
                
                // Children particles to the outgoing signal process particles. All outgoing particles to the signal
                // vertex share the same children, so the index used to access the vector is irrelevant, here being "0".
                vector<ConstGenParticlePtr> signal_children = (vertex -> particles_out())[0] -> children();
                
                for (auto & child : signal_children) {
                    if (child -> end_vertex() -> status() == 5) taken_ids.push_back(child -> id());
                }
                continue;
                
            }
            
            else if (v_stat == hcol_status){
                hcol_counter += 1;
                continue;
            }
            
            else if (v_stat == frag_status){
                for(auto & prtcl : vertex -> particles_in()){
                    
                    Particle p = Particle(prtcl);  // Extract particle as a Rivet Particle object.
                    
                    double pT = p.momentum().pT();  // Transverse momentum.
                    double deta = p.eta() - Zeta;  // Difference in pseudorapidity.
                    double dphi = p.phi() - Zphi;  // Difference in azimuthal angles.
                    
                    // Set dphi between (-pi, pi).
                    for(; std::fabs(dphi) > M_PI; dphi += (dphi > 0. ? -2.*M_PI : 2.*M_PI) );
                    
                    // Towards region
                    if(std::fabs(dphi) < M_PI/3.) pTsumTow += pT;
                    
                    // Transverse region
                    else if( std::fabs(dphi) < 2.*M_PI/3. ){
                        pTsumTrans += pT;
                        if(dphi > 0.) pTsumRight += pT;
                        else pTsumLeft += pT;
                    }
                    
                    // Away region
                    else pTsumAway += pT;

                    // Check if particle comes from the signal process vertex.
                    bool signal = false;
                    if(std::find(taken_ids.begin(), taken_ids.end(), prtcl -> id()) != taken_ids.end()) signal=true;
                    
                    // Store the values for writing them afterwards.
                    deta_values.push_back(deta);
                    dphi_values.push_back(dphi);
                    signals.push_back(signal);
                    pT_values.push_back(pT);
                    pid_values.push_back(p.pid());
                    
                }
                
                break;  // Here we assume that nothing of interest is after the Fragmentation vertex.
                
            }
        }
        
        // TransMAX, TransMIN regions
        pTsumTransmax = pTsumRight;
        pTsumTransmin = pTsumLeft;
        
        if (pTsumLeft > pTsumRight) {
            pTsumTransmax = pTsumLeft;
            pTsumTransmin = pTsumRight;
        }
        
        // Write down the collisions
        outfile << std::setprecision(my_width2 - 3);
        outfile << std::setw(my_width1) << std::left << "Hard_colls" << ", " << std::setw(my_width2) << std::right << hcol_counter << "\n";
        
        if (deta_values.size() == 0) return;  // Nothing more to be done if fragmentation vertex is empty.
        
        outfile << std::setprecision(my_width2 - 3);
        outfile << std::setw(my_width1) << std::left << "pTsumTow" << ", " << std::setw(my_width2) << std::right << pTsumTow << "\n";
        outfile << std::setw(my_width1) << std::left << "pTsumTra" << ", " << std::setw(my_width2) << std::right << pTsumTrans << "\n";
        outfile << std::setw(my_width1) << std::left << "pTsumAway" << ", " << std::setw(my_width2) << std::right << pTsumAway << "\n";
        outfile << std::setw(my_width1) << std::left << "pTsumTramin" << ", " << std::setw(my_width2) << std::right << pTsumTransmin << "\n";
        outfile << std::setw(my_width1) << std::left << "pTsumTramax" << ", " << std::setw(my_width2) << std::right << pTsumTransmax << "\n";
        write_vector(outfile, "pT", pT_values);
        write_vector(outfile, "dEta", deta_values);
        write_vector(outfile, "dPhi", dphi_values);
        write_vector(outfile, "Signal", signals);
        write_vector(outfile, "PIDs", pid_values);
        
    }
    // -----------------------------------------------------------------------------------------------------------------
    // End of "classify_particles" function.
    // -----------------------------------------------------------------------------------------------------------------
    
    
    // -----------------------------------------------------------------------------------------------------------------
    // Copy of ATLAS_2019_I1736531 analysis of the final state.
    // -----------------------------------------------------------------------------------------------------------------
    void analyze_final_state(const Event& event, ofstream& outfile, double Zphi, double Zpt){
        double pTsumTow(0.0), pTsumTrans(0.0), pTsumTransmin(0.0), pTsumTransmax(0.0);
        double pTsumAway(0.0), pTsumLeft(0.0), pTsumRight(0.0);
        
        const Cut& pcut = ( (Cuts::abspid != PID::SIGMAMINUS) && (Cuts::abspid != PID::SIGMAPLUS) &&
                            (Cuts::abspid != PID::XIMINUS)    && (Cuts::abspid != PID::OMEGAMINUS) );

        Particles particles = apply<ChargedFinalState>(event, "cfs").particlesByPt(Cuts::pT > 0.5*GeV && Cuts::abseta <2.5 && pcut);
        
        // Loop over charged particles
        
        for(const Particle& p : particles) {
            double dphi = p.momentum().phi() - Zphi;
            double pT   = p.momentum().pT();
            for(; std::fabs(dphi) > M_PI; dphi += (dphi > 0. ? -2.*M_PI : 2.*M_PI) );
            
            // Towards region
            if( std::fabs(dphi) < M_PI/3.) pTsumTow += pT;
            
            // Transverse region
            else if( std::fabs(dphi) < 2.*M_PI/3. ) {
                pTsumTrans += pT;
                if(dphi > 0.) pTsumRight += pT;
                else pTsumLeft += pT;
            }
            
            // Away region
            else pTsumAway += pT;
        }
    
        // TransMAX, TransMIN regions
        if (pTsumLeft > pTsumRight) {
            pTsumTransmax = pTsumLeft;
            pTsumTransmin = pTsumRight;
        }

        else {
            pTsumTransmax = pTsumRight;
            pTsumTransmin = pTsumLeft;
        }
        
        outfile << std::setprecision(my_width2-3);
        outfile << std::setw(my_width1) << std::left << "f_pTsumTow" << ", " << std::setw(my_width2) << std::right << pTsumTow/area << "\n";
        outfile << std::setw(my_width1) << std::left << "f_pTsumTra" << ", " << std::setw(my_width2) << std::right << pTsumTrans/area << "\n";
        outfile << std::setw(my_width1) << std::left << "f_pTsumAway" << ", " << std::setw(my_width2) << std::right << pTsumAway/area << "\n";
        outfile << std::setw(my_width1) << std::left << "f_pTsumTramin" << ", " << std::setw(my_width2) << std::right << pTsumTransmin/(0.5*area) << "\n";
        outfile << std::setw(my_width1) << std::left << "f_pTsumTramax" << ", " << std::setw(my_width2) << std::right << pTsumTransmax/(0.5*area) << "\n";
        
    }
    // -----------------------------------------------------------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------
    
    
    // -----------------------------------------------------------------------------------------------------------------
    // Global variables.
    // -----------------------------------------------------------------------------------------------------------------
    std::vector<int> Zpt_bins{0, 10, 20, 40, 60, 80, 120, 200};
    std::vector<string> data_paths;  // Store the paths to the directories of each Zpt range.
    const double area = 5.*2./3.*M_PI;
    int my_width1 = 15;  // Width for headers in output file.
    int my_width2 = 13;  // Width for values in output file.
    // -----------------------------------------------------------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------
    
    
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(FRAGMENTATION_CLASSIFICATION);
    
    void init() {
        
        // Get options from the new option system
        PdgId flav = (getOption("LMODE") == "EL")? PID::ELECTRON : PID::MUON;
        
        //Projections
        FinalState fs;
        ZFinder zfinder(
                        fs,  // Final state.
                        Cuts::abseta<2.4 && Cuts::pT>25.0*GeV,  //
                        flav,  // PdgId.
                        66*GeV,  // Minimum dilepton mass.
                        116*GeV,  // Maximum dilepton mass.
                        0.1,  // Maximum dR of photons around leptons to take into account for Z reconstruction.
                        ZFinder::ClusterPhotons::NODECAY  // Whether such photons are supposed to be clustered to the lepton objects and thus Z mom.
                        );
        declare(zfinder, "ZFinder");
        ChargedFinalState cfs(zfinder.remainingFinalState() );
        declare(cfs, "cfs");
        
        // -----------------------------------------------------------------------------------------------------------------
        // Create primary and secondary directories:
        // For each Zpt range, there is an associated secondary directory in which the particles files are saved.
        // The primary directory contains all the second directories.
        // -----------------------------------------------------------------------------------------------------------------
        const char* prim_dir = "./Fragmentation_classification_1M/";  // Primary directory for the output files.
        string path1 = prim_dir;

        // Create the primary directory if it doesn't exist.
        struct stat sb;  // Structure which would store the metadata
        if (!(stat(prim_dir, &sb) == 0)){
            std::__fs::filesystem::create_directory(prim_dir);
        }
        
        // -----------------------------------------------------------------------------------------------------------------
        // Create the secondary directories if they don't exist.
        for (int i = 0; i < Zpt_bins.size() - 1; i++) {
            
            string path2 = path1 + to_string(i+1) + "_" + to_string(Zpt_bins[i]) + " < Zpt < " + to_string(Zpt_bins[i + 1]);
            
            data_paths.push_back(path2);
            
            if (!(stat((path2).c_str(), &sb) == 0)){
                std::__fs::filesystem::create_directory(path2.c_str());
            }
        }
        
        string path2 = path1 + to_string(Zpt_bins.size()) + "_" + "Zpt >" + to_string(Zpt_bins[Zpt_bins.size()-1]);
        
        data_paths.push_back(path2);
        
        if (!(stat((path2).c_str(), &sb) == 0)){
            std::__fs::filesystem::create_directory(path2.c_str());
        }
        // -----------------------------------------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------------------------

    }
    
    
    // Perform the per-event analysis
    void analyze(const Event& event) {
        
        const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
        
        if (zfinder.bosons().size() != 1) vetoEvent;
        double  Zpt   = zfinder.bosons()[0].momentum().pT()/GeV;
        double  Zphi  = zfinder.bosons()[0].momentum().phi();
        double  Zeta  = zfinder.bosons()[0].momentum().eta();
        
        // -----------------------------------------------------------------------------------------------------------------
        // Determine Zpt region.
        // -----------------------------------------------------------------------------------------------------------------
        int i_bin(Zpt_bins.size() - 1);
                
        for (int i = 1; i < Zpt_bins.size(); i++) {
            if (Zpt < Zpt_bins[i]) {
                i_bin = i - 1;
                break;
            }
        }
        // -----------------------------------------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------------------------
        
        
        // -----------------------------------------------------------------------------------------------------------------
        // Create the output file in the corresponding directory.
        // -----------------------------------------------------------------------------------------------------------------
        string path = data_paths[i_bin];
        
        // Initialize counter and, if a file with the given value already exist,
        // increase the counter until this is no longer the case.
        int counter = 0;
        while (exists(path+"/particles_"+to_string(counter)+".csv")){
            counter++;
        }
        
        ofstream particle_file (path+"/particles_"+to_string(counter)+".csv");  // Open file;
        
        
        particle_file << std::setprecision(my_width2 - 3);
        particle_file << std::setw(my_width1) << std::left << "ZpT (GeV)" << ", " << std::right << std::setw(my_width2) << Zpt << "\n";
        particle_file << std::setw(my_width1) << std::left << "Z_eta" << ", " << std::setw(my_width2) << std::right << Zeta << "\n";
        
        // Perform ATLAS analysis.
        analyze_final_state(event, particle_file, Zphi, Zpt);
        
        // Classification of the incoming particles to fragmentation vertex.
        classify_particles(event, particle_file, Zphi, Zeta, Zpt);
        
        particle_file.close();
        
    }
    
    
};

// This global object acts as a hook for the plugin system
RIVET_DECLARE_PLUGIN(FRAGMENTATION_CLASSIFICATION);

}

