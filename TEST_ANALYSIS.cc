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


#include "Rivet/Event.hh"
#include "Rivet/Particle.hh"

namespace Rivet {

/// @brief Underlying event in Z events
class TEST_ANALYSIS : public Analysis {
public:
    
    inline bool exists (const std::string& filename) {
      struct stat buffer;
      return (stat (filename.c_str(), &buffer) == 0);
    }
    
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TEST_ANALYSIS);
    
    void init() {
        
        // Get options from the new option system
        PdgId flav = (getOption("LMODE") == "EL")? PID::ELECTRON : PID::MUON;
        
        //Projections
        FinalState fs;
        ZFinder zfinder(fs, Cuts::abseta<2.4 && Cuts::pT>25.0*GeV, flav, 66*GeV, 116*GeV, 0.1, ZFinder::ClusterPhotons::NODECAY);
        declare(zfinder, "ZFinder");
        ChargedFinalState cfs(zfinder.remainingFinalState() );
        declare(cfs, "cfs");
    }
    
    
//    // Victor's function to get the first vertex for a given signal value
//    shared_ptr<const HepMC3::GenVertex> get_single_vertex(int my_stat, const std::vector<ConstGenVertexPtr> the_vertices)
//    {
//        const auto & iterator = std::find_if(the_vertices.begin(), the_vertices.end(), [my_stat](const std::shared_ptr<const HepMC3::GenVertex> v)
//                                                                                  {
//                                                                                     return (v -> status()) == my_stat;
//                                                                                  });
//        return *iterator;
//    }
    
    
    // Perform the per-event analysis
    void analyze(const Event& event) {
        
        const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
        
        cout << string(120, '=') << "\n" << endl;
        cout << "MY OUTPUTS GO BELOW" << std::endl;
        cout << "\n" << endl;
        
        // ------------------------------------------------------------------------------------------------------------------------------------------
        // Basic tests on handling HepMC3. Constant version.
        // ------------------------------------------------------------------------------------------------------------------------------------------
        
        // The variable "event" defined at the beginning of "analyze" is a Rivet variable.
        // It can be un-wrapped to a HepMC3 object as follows:
        
        const GenEvent evt = *event.genEvent();  // Constant version.
        
        // The properties of the new variable "evt" can be read in https://gitlab.cern.ch/hepmc/HepMC3/-/blob/master/include/HepMC3/GenEvent.h
        // Of interest for me, is to access the vertices of the event.
        
        const std::vector<ConstGenVertexPtr> vertices = evt.vertices();
        
        // The defined "vertices" variable is a vector and its elements can be accessed just like 'vertices[1]'.
        
//        RivetHepMC::Print::line(vertices[1]);  // Print in a nice format the first vertex in vertices.
        
        // One can also extract one of the vertices via its index as follows:
        
//        const std::shared_ptr<const HepMC3::GenVertex> &first_vertex = vertices[1];
        
//        RivetHepMC::Print::line(first_vertex);  // Print in a nice format the first vertex in vertices.
        
        // A simple for-loop to print every vertex in 'vertices' can be written as follows.
//        for(auto & vertex : vertices){
//            RivetHepMC::Print::line(vertex);
//        }
        
        // Each vertex in "vertices" has an attribute "status", which is an integer representing some specific process in that vertex.
        // Using Sherpa, some status numbers of interest are: (1, Signal Process), (5, Fragmentation) and (2, Hard Collision).
        // One can then extract a (single) vertex by its status as follows:
        
        //   1. Define the number of status of interest.
        int get_status = 1;
        
        //   2. Define the vertex variable in which to save the extracted vertex.
        
        shared_ptr<const HepMC3::GenVertex> my_vertex;
        
        //   3. Copy the vertex from "vertices".
        
        for (auto & vertex : vertices) {
            if((vertex -> status()) == get_status){
                my_vertex = vertex;
                break;
            }
        }
        
//        RivetHepMC::Print::line(my_vertex);  // Print in a nice format the resulting vertex.
        
        
        // Any vertex has the properties "particles_in" and "particles_out". These are vectors
        // with every entry corresponding to a particle. They can be accessed as follows:
        
        const std::vector<ConstGenParticlePtr> my_particles = my_vertex -> particles_out();
        
        // A simple for-loop to print every particle in 'my_particles'.
//        for(auto & particle : my_particles){
//            RivetHepMC::Print::line(particle);
//        }
        
        // Once the incoming/outgoing particles of the chosen vertex are taken, one can define another
        // vector of particles to be filled from the first vector, according to certain criteria.
        // Here, for example, from the childrens of the particles of the signal process, those which
        // end in the fragmentation are taken.
        
        //   1. Define the children particles.
        // Here we access the first particle from "my_particles" because, apparently, any
        // particle here has the same children, so which index is used makes no difference.

        vector<ConstGenParticlePtr> children = my_particles[1] -> children();
        
        //   2. Define the empty vector for the chosen particles and proceed to fill it
        //   with the children particles whose end vertex is the fragmentation vertex.
        
        vector<shared_ptr<const HepMC3::GenParticle>> chosen_particles;
        
        std::copy_if(children.begin(), children.end(), std::back_inserter(chosen_particles), [](std::shared_ptr<const HepMC3::GenParticle> child)
                     {
            return (child -> end_vertex()) -> status() == 5;
            
        });
        
        // Create a vector with the ids of the chosen particles for future distinction.
        std::vector<int> taken_ids;
        std::transform(chosen_particles.begin(), chosen_particles.end(), std::back_inserter(taken_ids), [](std::shared_ptr<const HepMC3::GenParticle> chosen_p)
                       {
            return chosen_p -> id();

        });
        
        // Get the Fragmentation vertex and its incoming particles.
        shared_ptr<const HepMC3::GenVertex> fragmentation_vertex;

        for (auto & vertex : vertices) {
            if((vertex -> status()) == 5){
                fragmentation_vertex = vertex;
                break;
            }
        }
        
        const vector<shared_ptr<const HepMC3::GenParticle>> frag_particles = fragmentation_vertex -> particles_in();
        
        
        // Get those particles incoming to the fragmentation vertex that are not children of the outgoing signal process particles.

        vector<shared_ptr<const HepMC3::GenParticle>> not_chosen;

        std::copy_if(frag_particles.begin(), frag_particles.end(), std::back_inserter(not_chosen), [taken_ids](std::shared_ptr<const HepMC3::GenParticle> frag_p)
                     {
            return !(std::find(taken_ids.begin(), taken_ids.end(), frag_p -> id()) != taken_ids.end());

        });
        
        cout << "\n" << endl;
        cout << "All fragmentation incoming particles" << endl;
        cout << "\n" << endl;


        for(auto & particle : frag_particles){
            RivetHepMC::Print::line(particle);

        }

        cout << "\n" << endl;
        cout << "Those children to the signal process" << endl;
        cout << "\n" << endl;


        for(auto & particle : chosen_particles){
            RivetHepMC::Print::line(particle);

        }

        cout << "\n" << endl;
        cout << "Those due to other processes" << endl;
        cout << "\n" << endl;


        for(auto & particle : not_chosen){
            RivetHepMC::Print::line(particle);
//            cout << Particle(particle).phi() << endl;
        }
        
        
        // ------------------------------------------------------------------------------------------------------------------------------------------
        // ------------------------------------------------------------------------------------------------------------------------------------------
        
        
        
        
        
//        // ------------------------------------------------------------------------------------------------------------------------------------------
//        // Non-constant version
//        // ------------------------------------------------------------------------------------------------------------------------------------------
//
//        // Get event as HepMC3 variable.
//        GenEvent evt = *event.genEvent();
//
//        // Extract vertices.
//        std::vector<HepMC3::GenVertexPtr> vertices = evt.vertices();
//
//        // Extract the signal process vertex, which has status = 1.
//        int get_status = 1;
//
//        shared_ptr<HepMC3::GenVertex> signal_vertex;
//
//        for (auto & vertex : vertices) {
//            if((vertex -> status()) == get_status){
//                signal_vertex = vertex;
//                break;
//            }
//        }
//
//        // Get the outgoing particles from the signal process vertex.
//        std::vector<HepMC3::GenParticlePtr> my_particles = signal_vertex -> particles_out();
//
//        // Get the children of the outgoing particles from the signal process vertex.
//        std::vector<HepMC3::GenParticlePtr> children = my_particles[1] -> children();
//
//        // Chose the children that end in the Fragmentation vertex.
//        std::vector<HepMC3::GenParticlePtr> chosen_particles;
//
//        std::copy_if(children.begin(), children.end(), std::back_inserter(chosen_particles), [](std::shared_ptr<HepMC3::GenParticle> child)
//                     {
//            return (child -> end_vertex()) -> status() == 5;
//
//        });
//
//
//        // Create a vector with the ids of the chosen particles.
//        std::vector<int> taken_ids;
//        std::transform(chosen_particles.begin(), chosen_particles.end(), std::back_inserter(taken_ids), [](std::shared_ptr<HepMC3::GenParticle> chosen_p)
//                       {
//            return chosen_p -> id();
//
//        });
//
//
//        // Get the Fragmentation vertex and its incoming particles.
//        shared_ptr<HepMC3::GenVertex> fragmentation_vertex;
//
//        for (auto & vertex : vertices) {
//            if((vertex -> status()) == 5){
//                fragmentation_vertex = vertex;
//                break;
//            }
//        }
//
//        std::vector<HepMC3::GenParticlePtr> frag_particles = fragmentation_vertex -> particles_in();
//
//
//        // Get those particles incoming to the fragmentation vertex that are not children of the outgoing signal process particles.
//
//        vector<shared_ptr<HepMC3::GenParticle>> not_chosen;
//
//        std::copy_if(frag_particles.begin(), frag_particles.end(), std::back_inserter(not_chosen), [taken_ids](std::shared_ptr<HepMC3::GenParticle> frag_p)
//                     {
//            return !(std::find(taken_ids.begin(), taken_ids.end(), frag_p -> id()) != taken_ids.end());
//
//        });
//
//
//        cout << "\n" << endl;
//        cout << "All fragmentation incoming particles" << endl;
//        cout << "\n" << endl;
//
//
//        for(auto & particle : frag_particles){
//            RivetHepMC::Print::line(particle);
//
//        }
//
//        cout << "\n" << endl;
//        cout << "Those children to the signal process" << endl;
//        cout << "\n" << endl;
//
//
//        for(auto & particle : chosen_particles){
//            RivetHepMC::Print::line(particle);
//
//        }
//
//        cout << "\n" << endl;
//        cout << "Those due to other processes" << endl;
//        cout << "\n" << endl;
//
//
//        for(auto & particle : not_chosen){
//            RivetHepMC::Print::line(particle);
//
//        }
//
//
//        // ------------------------------------------------------------------------------------------------------------------------------------------
//        // ------------------------------------------------------------------------------------------------------------------------------------------
        
        
        
        cout << "\n" << endl;
        cout << "END OF MY OUTPUTS" << std::endl;
        cout << string(120, '=') << "\n" << endl;
        
//        if (zfinder.bosons().size() != 1) vetoEvent;
//        double  Zphi  = zfinder.bosons()[0].momentum().phi();
//        double  Zeta  = zfinder.bosons()[0].momentum().eta();

    
//        std::vector<double> leftpt;
//        std::vector<double> rightpt;
//
//        const Cut& pcut = ( (Cuts::abspid != PID::SIGMAMINUS) && (Cuts::abspid != PID::SIGMAPLUS) &&
//                           (Cuts::abspid != PID::XIMINUS)    && (Cuts::abspid != PID::OMEGAMINUS) );
//
//        Particles particles = apply<ChargedFinalState>(event, "cfs").particlesByPt(Cuts::pT > 0.5*GeV && Cuts::abseta <2.5 && pcut);
//
//
//        //Calculate thrust
//        vector<Vector3> momenta;
//        for(const Particle& p : particles) {
//            Vector3 mom = p.momentum().vector3();
//            mom.setZ(0.0);
//            momenta.push_back(mom);
//        }
//
//        if (momenta.size() == 2) {
//            momenta.push_back(Vector3(1e-10*MeV, 0., 0.));
//        }
//
//        Thrust thrustC;
//        thrustC.calc(momenta);
//
//
//        // Start writing the path of the file in which the data will be written. Folder must exist beforehand.
//
//        string path = "./new_EventDisplay/";
//
//        // Initialize the counter and, if a file with the given value already exist, increase the counter untill this is no longer the case.
//        int counter = 0;
//        while (exists(path+"particles("+to_string(counter)+").csv")){
//            counter++;
//        }
//        ofstream particle_file (path+"particles("+to_string(counter)+").csv");
//
//        particle_file << "dEta\tdPhi" << endl;  // Header of the file.
//
//        // Loop over charged particles
//        for(const Particle& p : particles) {
//            double dphi = p.momentum().phi() - Zphi;
//            double deta = p.momentum().eta() - Zeta;
//            for(; std::fabs(dphi) > M_PI; dphi += (dphi > 0. ? -2.*M_PI : 2.*M_PI) );
//
//            particle_file << to_string(deta) << "\t" << to_string(dphi) << std::endl;
//        }
//
//        particle_file.close();
        
    }
    
    
};

// This global object acts as a hook for the plugin system
RIVET_DECLARE_PLUGIN(TEST_ANALYSIS);

}
