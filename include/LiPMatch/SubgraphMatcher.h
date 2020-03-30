#ifndef __SUBGRAPHMATCHER_H
#define __SUBGRAPHMATCHER_H



#include <Subgraph.h>
#include <Plane.h>
#include <Pole.h>
#include <Vehicle.h>

namespace LiPMatch_ns {


  class SubgraphMatcher
  {
   public:

    SubgraphMatcher();

    bool evalUnaryConstraintsPlane(Plane &plane1, Plane &plane2);

    bool evalUnaryConstraintsVehicle(Vehicle &vehicle1, Vehicle &vehicle2);

    bool evalUnaryConstraintsPole(Pole &pole1, Pole &pole2);

    bool evalBinaryConstraints(Plane &plane1, Plane &plane2, Plane &planeA, Plane &planeB);

    bool evalBinaryConstraintsVehicle(Vehicle &Ref, Vehicle &neigRef, Vehicle &Check, Vehicle &neigCheck);

    bool evalBinaryConstraintsPole(Pole &Ref, Pole &neigRef, Pole &Check, Pole &neigCheck);

    std::vector<std::map<unsigned,unsigned> > alreadyExplored;

    void exploreSubgraphTreeR(std::set<unsigned> &evalRef, std::set<unsigned> &evalCheck, std::map<unsigned, unsigned> &matched);

    void exploreSubgraphTreeRVehicle(std::set<unsigned> &sourceVehicles, std::set<unsigned> &targetVehicles, std::map<unsigned, unsigned> &matched);

    void exploreSubgraphTreeRPole(std::set<unsigned> &sourcePoles, std::set<unsigned> &targetPoles, std::map<unsigned, unsigned> &matched);



    std::map<unsigned,unsigned> compareSubgraphs(Subgraph &subgraphSource, Subgraph &subgraphTarget, int& unaryCount);


    std::map<unsigned,unsigned> compareSubgraphsVehiclePlaneRef(Subgraph &subgraphSource, Subgraph &subgraphTarget, int& unaryCount,
                                                                std::vector<Eigen::Vector3f>& kvc, std::vector<Eigen::Vector3f>& lvc,
                                                                std::vector<Eigen::Vector3f>& kvn, std::vector<Eigen::Vector3f>& lvn);

    std::map<unsigned,unsigned> compareSubgraphsPolePlaneRef(Subgraph &subgraphSource, Subgraph &subgraphTarget, int& unaryCount,
                                                             std::vector<Eigen::Vector3f>& kvc, std::vector<Eigen::Vector3f>& lvc,
                                                             std::vector<Eigen::Vector3f>& kvn, std::vector<Eigen::Vector3f>& lvn);


    Subgraph *subgraphSrc;

    Subgraph *subgraphTrg;

    float wdif_height;
    float wdif_height2;
    float wdif_normal;
    float wrel_dist_centers;

    float wal;
    float wea;


//    float height_threshold = 1.75;
//    float angle_threshold = 1.75;
//    float dist_threshold = 1.75;
//
////      matcher.height_threshold = 1.25;
////      matcher.angle_threshold = 2.0;
////      matcher.dist_threshold = 1.5;
//
//    float area_threshold = 1.65;
//    float elongation_threshold = 1.65;

      float radios = 1.72;

      float height_threshold = radios;
      float angle_threshold = radios;
      float dist_threshold = radios;
      float area_threshold = radios;
      float elongation_threshold = radios;


//        float height_threshold = 1.7;
//        float angle_threshold = 1.7;
//        float dist_threshold = 1.7;
//
////      matcher.height_threshold = 1.25;
////      matcher.angle_threshold = 2.0;
////      matcher.dist_threshold = 1.5;
//
//        float area_threshold = 1.7;
//        float elongation_threshold = 1.7;


  private:

    std::map<unsigned, unsigned> winnerMatch;

    std::map<unsigned, unsigned> winnerMatchVehicle;

    std::map<unsigned, unsigned> winnerMatchPole;

    std::vector<std::vector<int8_t> > hashUnaryConstraints;

    std::vector<Eigen::Vector3f> c_kvc;
    std::vector<Eigen::Vector3f> c_lvc;
    std::vector<Eigen::Vector3f> c_kvn;
    std::vector<Eigen::Vector3f> c_lvn;

    std::vector<float> vdif_height;
    std::vector<float> vdif_height2;
    std::vector<float> vdif_normal;
    std::vector<float> vrel_dist_centers;

    std::vector<float> vah;
    std::vector<float> vea;

  };

 } // End of namespaces

#endif
