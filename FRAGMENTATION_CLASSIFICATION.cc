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
    
    // -----------------------------------------------------------------------------------------------------------------
    // Function to classify the particles incoming to the Fragmentation vertex according to
    // whether or not they are children to the outgoing particles of the signal process.
    // -----------------------------------------------------------------------------------------------------------------
    void classify_particles(const Event& event, ofstream& outfile, double area,
                             double Zphi, double Zeta, double Zpt, double Zthe, int h_length) {
        int sign_status = 1;  // Signal process vertex has status = 1.
        int hcol_status = 2;  // Hard collision vertices have status = 2.
        int frag_status = 5;  // Fragmentation vertex has status = 5.
        int hcol_counter = 0;  // Hard collisions counter.
        
        std::vector<int> taken_ids;  // store ids of the children to the outgoing particles of the signal vertex.
        
        std::vector<float> deta_values, dphi_values, dthe_values;
        std::vector<int> signals;
        
        double pTsumTow(0.0), pTsumTrans(0.0), pTsumAway(0.0), pTsum(0.0);
        double pTsumTransmin(0.0), pTsumTransmax(0.0), pTsumLeft(0.0), pTsumRight(0.0);
        
        ///@todo: Implement cutting Cuts::pT > 0.5*GeV && Cuts::abseta <2.5 && pcut
        
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
                // Testing testing
                for(auto & p : vertex -> particles_in()){
                    double deta = Particle(p).eta() - Zeta;  // Difference in pseudorapidity.
                    double dphi = Particle(p).phi() - Zphi;  // Difference in azimuthal angles.
                    double dthe = Particle(p).theta() - Zthe;  // Difference in polar angles.
                    double pT = Particle(p).momentum().pT();  // Transverse momentum.
                    
                    pTsum += pT;
                    
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

                    bool signal = false;
                    // If true, particle comes from the signal process.
                    if(std::find(taken_ids.begin(), taken_ids.end(), p -> id()) != taken_ids.end()) signal=true;
                    
                    // Store the values for writing them afterwards.
                    deta_values.push_back(deta);
                    dphi_values.push_back(dphi);
                    dthe_values.push_back(dthe);
                    signals.push_back(signal);
                    
                }
                
                break;  // Here we assume that nothing of interest is after the Fragmentation vertex.
                
            }
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
        
        for (int i = 0; i < deta_values.size() - 1; i++){
            outfile << std::setprecision(4) << std::fixed;
            outfile << std::setw(8) << deta_values[i] << "\t";
            outfile << std::setw(8) << dphi_values[i] << "\t";
            outfile << std::setw(8) << dthe_values[i] << "\t";
            
            if(i != 0){
                outfile << std::setw(7) << signals[i] << "\n";
                continue;
            }
            
            outfile << std::setw(7) << signals[i] << "\t";
            outfile << std::setw(8) << Zpt << "\t";
            outfile << std::setw(8) << hcol_counter << "\t";
            outfile << std::setw(8) << pTsum << "\t";
            outfile << std::setw(8) << pTsumTow << "\t";
            outfile << std::setw(10) << pTsumTrans << "\t";
            outfile << std::setw(9) << pTsumAway << "\t";
            outfile << std::setw(13) << pTsumTransmin << "\t";
            outfile << std::setw(13) << pTsumTransmax << "\t";
            outfile << std::setw(9) << area << "\t";
            
            analyze_final_state(event, outfile, area, Zphi, Zpt);
//            outfile << std::setw(9) << to_string(area) << '\n';
        }
        
    }
    // -----------------------------------------------------------------------------------------------------------------
    // End of "classify_particles" function.
    // -----------------------------------------------------------------------------------------------------------------
    
    
    // -----------------------------------------------------------------------------------------------------------------
    // Copy of ATLAS_2019_I1736531 analysis of the final state.
    // -----------------------------------------------------------------------------------------------------------------
    void analyze_final_state(const Event& event, ofstream& outfile, double area, double Zphi, double Zpt){
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
        
        outfile << std::setprecision(4) << std::fixed;
        outfile << std::setw(10) << pTsumTow/area << "\t";
        outfile << std::setw(12) << pTsumTrans/area << "\t";
        outfile << std::setw(11) << pTsumAway/area << "\t";
        outfile << std::setw(15) << pTsumTransmin/(0.5*area) << "\t";
        outfile << std::setw(15) << pTsumTransmax/(0.5*area) << '\n';
        
    }
    // -----------------------------------------------------------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------
    
    
    // -----------------------------------------------------------------------------------------------------------------
    // Global variables.
    // -----------------------------------------------------------------------------------------------------------------
    std::vector<int> Zpt_bins{0, 10, 20, 40, 60, 80, 120, 200};
    int files_per_value = 10;  // Maximum amount of files to generate per bin of Zpt.
    std::vector<string> data_paths;  // Store the paths to the directories of each Zpt range.
    std::vector<int> data_counter = std::vector<int>(Zpt_bins.size());  // Counts the files generated per Zpt entry.
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
        const char* prim_dir = "./Fragmentation_classification2/";  // Primary directory for the output files.
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
        
        const double area = 5.*2./3.*M_PI;
        const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
        
        if (zfinder.bosons().size() != 1) vetoEvent;
        double  Zpt   = zfinder.bosons()[0].momentum().pT()/GeV;
        double  Zphi  = zfinder.bosons()[0].momentum().phi();
        double  Zeta  = zfinder.bosons()[0].momentum().eta();
        double  Zthe  = zfinder.bosons()[0].momentum().theta();
        
        ///@todo: Implement a function that creates statistics such that one can extract representatives from the sample.
        
        // -----------------------------------------------------------------------------------------------------------------
        // Determine Zpt region.
        // -----------------------------------------------------------------------------------------------------------------
        int i_bin(0);
        
        for (int i = 0; i < Zpt_bins.size() - 1; i++) {
            if (inRange(Zpt,Zpt_bins[i], Zpt_bins[i+1])){
                i_bin = i;
                break;
            }
        }
        
        if (Zpt > Zpt_bins[Zpt_bins.size() - 1]) i_bin = Zpt_bins.size() - 1;
        // -----------------------------------------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------------------------
        
        // If the corresponding region already has enough files, veto the event.
//        if (data_counter[i_bin] > files_per_value) vetoEvent;
        
//        data_counter[i_bin] += 1;
        
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
        
        // Header of the file
        particle_file << std::setw(8) << "dEta" << "\t";
        particle_file << std::setw(8) << "dPhi" << "\t";
        particle_file << std::setw(8) << "dTheta" << "\t";
        particle_file << std::setw(7) << "S_child" << "\t";
        particle_file << std::setw(8) << "ZpT" << "\t";
        particle_file << std::setw(8) << "Hard_col" << "\t";
        particle_file << std::setw(8) << "pTsum" << "\t";
        particle_file << std::setw(8) << "pTsumTow" << "\t";
        particle_file << std::setw(10) << "pTsumTrans" << "\t";
        particle_file << std::setw(9) << "pTsumAway" << "\t";
        particle_file << std::setw(13) << "pTsumTransmin" << "\t";
        particle_file << std::setw(13) << "pTsumTransmax" << "\t";
        particle_file << std::setw(9) << "dEta*dPhi" << "\t";
        particle_file << std::setw(10) << "f_pTsumTow" << "\t";
        particle_file << std::setw(12) << "f_pTsumTrans" << "\t";
        particle_file << std::setw(11) << "f_pTsumAway" << "\t";
        particle_file << std::setw(15) << "f_pTsumTransmin" << "\t";
        particle_file << std::setw(15) << "f_pTsumTransmax" << '\n';
        
        int h_length = 8 + 8 + 8 + 7 + 8 + 8 + 8 + 8 + 10 + 9 + 13 + 13 + 9 + 10 + 12 + 11 + 15 + 15;
        
        classify_particles(event, particle_file, area, Zphi, Zeta, Zpt, Zthe, h_length);
        
        particle_file.close();
        
    }
    
    
};

// This global object acts as a hook for the plugin system
RIVET_DECLARE_PLUGIN(FRAGMENTATION_CLASSIFICATION);

}

