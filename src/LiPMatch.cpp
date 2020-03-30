#include <LiPMatch.h>

using namespace std;
using namespace Eigen;
using namespace LiPMatch_ns;


LiPMatch::LiPMatch() : LiPMatch_stop(false), LiPMatch_finished(false)
{
    LiPMatch_hd = createThreadFromObjectMethod(this,&LiPMatch::run);

    keyframe_vec.clear();

    all_kfs_pose.points.clear();

    map_rfn.set_down_sample_resolution( 0.75 );
}



void LiPMatch::detectPlanesCloud( m_keyframe &c_keyframe, int keyFrameCount)
{
//    static std::ofstream unarycout("/home/jjwen/u.txt");
//    static std::ofstream binarycout("/home/jjwen/b.txt");

    static pcl::VoxelGrid<pcl::PointXYZI> grid_1;
    grid_1.setLeafSize(0.35,0.35,0.35);
//    grid_1.setInputCloud(c_keyframe.orilaserCloud);
//    grid_1.filter (*c_keyframe.orilaserCloud);

    grid_1.setInputCloud(c_keyframe.structurelaserCloud);
    grid_1.filter (*c_keyframe.structurelaserCloud);

    grid_1.setInputCloud(c_keyframe.vehiclelaserCloud);
    grid_1.filter (*c_keyframe.vehiclelaserCloud);

    grid_1.setInputCloud(c_keyframe.naturelaserCloud);
    grid_1.filter (*c_keyframe.naturelaserCloud);

    grid_1.setInputCloud(c_keyframe.objectlaserCloud);
    grid_1.filter (*c_keyframe.objectlaserCloud);

    static pcl::VoxelGrid<pcl::PointXYZI> grid_2;
    grid_2.setLeafSize(0.8,0.8,0.8);
    grid_2.setInputCloud(c_keyframe.surflaserCloud);
    grid_2.filter (*c_keyframe.surflaserCloud);

    static pcl::VoxelGrid<pcl::PointXYZI> grid_3;
    grid_3.setLeafSize(0.4,0.4,0.4);
    grid_3.setInputCloud(c_keyframe.linelaserCloud);
    grid_3.filter (*c_keyframe.linelaserCloud);

    grid_3.setInputCloud(c_keyframe.orilaserCloud);
    grid_3.filter (*c_keyframe.orilaserCloud);

    grid_3.setInputCloud(c_keyframe.g_laserCloud.makeShared());
    grid_3.filter (c_keyframe.g_laserCloud);

    double time_behinSeg = pcl::getTime();

    //////////////////////
    int colorid = 0;
    //聚类
    pcl::search::KdTree<pcl::PointXYZI>::Ptr treeVehicle (new pcl::search::KdTree<pcl::PointXYZI>);
    treeVehicle->setInputCloud (c_keyframe.vehiclelaserCloud);
    std::vector<pcl::PointIndices> cluster_indices_vehicle;
    pcl::EuclideanClusterExtraction<pcl::PointXYZI> ec_vehicle;
    ec_vehicle.setClusterTolerance (0.45);
//    ec.setMinClusterSize (100);
    ec_vehicle.setMinClusterSize (100);
    ec_vehicle.setMaxClusterSize (150000);
    ec_vehicle.setSearchMethod (treeVehicle);
    ec_vehicle.setInputCloud (c_keyframe.vehiclelaserCloud);
    ec_vehicle.extract (cluster_indices_vehicle);
    pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudNgAftEe_vehicles(new pcl::PointCloud<pcl::PointXYZI>());
    pcl::PCA< pcl::PointXYZI > pca;

    vector<Vehicle> detectedLocalVehicles;

    for (std::vector<pcl::PointIndices>::const_iterator it = cluster_indices_vehicle.begin (); it != cluster_indices_vehicle.end (); ++it)
    {
        pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudNgAftEe_vehicle(new pcl::PointCloud<pcl::PointXYZI>());
        for (std::vector<int>::const_iterator pit = it->indices.begin (); pit != it->indices.end (); ++pit)
        {
            c_keyframe.vehiclelaserCloud->points[*pit].intensity = colorid;
            laserCloudNgAftEe_vehicle->points.push_back (c_keyframe.vehiclelaserCloud->points[*pit]);
        }
        pca.setInputCloud(laserCloudNgAftEe_vehicle);
        Eigen::VectorXf eigenVal = pca.getEigenValues();

        if (eigenVal[2] / eigenVal[0] < 0.005)
            continue;

        *laserCloudNgAftEe_vehicles += *laserCloudNgAftEe_vehicle;

        colorid++;

//        std::cout<<eigenVal[0]<<" "<<eigenVal[1]<<" "<<eigenVal[2]<<std::endl;

        Vehicle vc;
        vc.VehiclePointCloudPtr = laserCloudNgAftEe_vehicle;
        vc.keyFrameId = keyFrameCount;
        vc.calcCenterAndElongation();

        detectedLocalVehicles.push_back(vc);

    }

    laserCloudOri_mp1 += *laserCloudNgAftEe_vehicles;

    observedVehicles.clear();
    vVehicles.clear();

    for (size_t i = 0; i < detectedLocalVehicles.size (); i++)
    {
        detectedLocalVehicles[i].id = vVehicles.size();

        // Update co-visibility graph
        for(set<unsigned>::iterator it = observedVehicles.begin(); it != observedVehicles.end(); it++)
        {
            detectedLocalVehicles[i].neighborVehicles[*it] = 1;
            vVehicles[*it].neighborVehicles[detectedLocalVehicles[i].id] = 1;
        }
        observedVehicles.insert(detectedLocalVehicles[i].id);
        vVehicles.push_back(detectedLocalVehicles[i]);
    }


    std::shared_ptr<Maps_keyframe<float>> mk = std::make_shared<Maps_keyframe<float>>();

    std::cout<<"vVehicles.size() "<<vVehicles.size()<<std::endl;

    mk->vVehicles = vVehicles;

    //聚类
    pcl::search::KdTree<pcl::PointXYZI>::Ptr treeNature (new pcl::search::KdTree<pcl::PointXYZI>);
    treeNature->setInputCloud (c_keyframe.naturelaserCloud);
    std::vector<pcl::PointIndices> cluster_indices_nature;
    pcl::EuclideanClusterExtraction<pcl::PointXYZI> ec_nature;
    ec_nature.setClusterTolerance (0.45);
//    ec.setMinClusterSize (100);
    ec_nature.setMinClusterSize (15);
    ec_nature.setMaxClusterSize (150000);
    ec_nature.setSearchMethod (treeNature);
    ec_nature.setInputCloud (c_keyframe.naturelaserCloud);
    ec_nature.extract (cluster_indices_nature);
    pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudNgAftEe_natures(new pcl::PointCloud<pcl::PointXYZI>());
    pcl::PCA< pcl::PointXYZI > pca1;

    vector<Pole> detectedLocalPoles;

    colorid = 0;
    for (std::vector<pcl::PointIndices>::const_iterator it = cluster_indices_nature.begin (); it != cluster_indices_nature.end (); ++it)
    {
        pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudNgAftEe_nature(new pcl::PointCloud<pcl::PointXYZI>());
        for (std::vector<int>::const_iterator pit = it->indices.begin (); pit != it->indices.end (); ++pit)
        {
            c_keyframe.naturelaserCloud->points[*pit].intensity = colorid;
            laserCloudNgAftEe_nature->points.push_back (c_keyframe.naturelaserCloud->points[*pit]);
        }
        pca1.setInputCloud(laserCloudNgAftEe_nature);
        Eigen::VectorXf eigenVal = pca1.getEigenValues();

        if (eigenVal[1] / eigenVal[0] > 0.17)
            continue;

        *laserCloudNgAftEe_natures += *laserCloudNgAftEe_nature;
        colorid++;

        Pole vc;
        vc.PolePointCloudPtr = laserCloudNgAftEe_nature;
        vc.keyFrameId = keyFrameCount;
        vc.calcCenterAndElongation();

        detectedLocalPoles.push_back(vc);

//        std::cout<<eigenVal[0]<<" "<<eigenVal[1]<<" "<<eigenVal[2]<<std::endl;

    }

//    laserCloudOri_mp1 += *laserCloudNgAftEe_natures;

    observedPoles.clear();
    vPoles.clear();

    for (size_t i = 0; i < detectedLocalPoles.size (); i++)
    {
        detectedLocalPoles[i].id = vPoles.size();

        // Update co-visibility graph
        for(set<unsigned>::iterator it = observedPoles.begin(); it != observedPoles.end(); it++)
        {
            detectedLocalPoles[i].neighborPoles[*it] = 1;
            vPoles[*it].neighborPoles[detectedLocalPoles[i].id] = 1;
        }
        observedPoles.insert(detectedLocalPoles[i].id);
        vPoles.push_back(detectedLocalPoles[i]);
    }

    std::cout<<"vPoles.size() "<<vPoles.size()<<std::endl;

    mk->vPoles = vPoles;


    static pcl::VoxelGrid<pcl::PointXYZI> grid_t;
    grid_t.setLeafSize(0.6,0.6,0.6);

    grid_t.setInputCloud(laserCloudOri_mp1.makeShared());
    grid_t.filter (laserCloudOri_mp1);



    /////////////////////



    //聚类
    pcl::search::KdTree<pcl::PointXYZI>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZI>);
    tree->setInputCloud (c_keyframe.structurelaserCloud);
    std::vector<pcl::PointIndices> cluster_indices;
    pcl::EuclideanClusterExtraction<pcl::PointXYZI> ec;
    ec.setClusterTolerance (0.8);
//    ec.setMinClusterSize (100);
    ec.setMinClusterSize (200);
    ec.setMaxClusterSize (150000);
    ec.setSearchMethod (tree);
    ec.setInputCloud (c_keyframe.structurelaserCloud);
    ec.extract (cluster_indices);
    pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudNgAftEe(new pcl::PointCloud<pcl::PointXYZI>());
    for (std::vector<pcl::PointIndices>::const_iterator it = cluster_indices.begin (); it != cluster_indices.end (); ++it)
    {
        for (std::vector<int>::const_iterator pit = it->indices.begin (); pit != it->indices.end (); ++pit)
        {
            laserCloudNgAftEe->points.push_back (c_keyframe.structurelaserCloud->points[*pit]);
        }
    }

//    pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudNgAftEe(new pcl::PointCloud<pcl::PointXYZI>());
//    laserCloudNgAftEe = c_keyframe.structurelaserCloud;

    //区域增长法提取激光点云中的平面
    pcl::search::Search<pcl::PointXYZI>::Ptr treeRg (new pcl::search::KdTree<pcl::PointXYZI>);
    pcl::PointCloud <pcl::Normal>::Ptr normals (new pcl::PointCloud <pcl::Normal>);
    pcl::NormalEstimation<pcl::PointXYZI, pcl::Normal> normal_estimator;
    normal_estimator.setSearchMethod (treeRg);
    normal_estimator.setInputCloud (laserCloudNgAftEe);
    normal_estimator.setKSearch (10);
    normal_estimator.compute (*normals);

    pcl::RegionGrowing<pcl::PointXYZI, pcl::Normal> reg;
//    reg.setMinClusterSize (100);
    reg.setMinClusterSize (100);
    reg.setMaxClusterSize (1000000);
    reg.setSearchMethod (treeRg);
    reg.setNumberOfNeighbours (20);
    reg.setInputCloud (laserCloudNgAftEe);
    reg.setInputNormals (normals);
    reg.setSmoothnessThreshold (14.0 / 180.0 * M_PI);
    reg.setCurvatureThreshold (10.0);

    std::vector <pcl::PointIndices> clusters;
    reg.extract (clusters);

    pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudNgAftReg(new pcl::PointCloud<pcl::PointXYZI>());
    //创建一个模型参数对象，用于记录结果
    pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
    //inliers表示误差能容忍的点 记录的是点云的序号
    pcl::PointIndices::Ptr inliers(new pcl::PointIndices);
    // 创建一个分割器
    pcl::SACSegmentation<pcl::PointXYZI> seg;
    // Optional，这个设置可以选定结果平面展示的点是分割掉的点还是分割剩下的点。
    seg.setOptimizeCoefficients(true);
    // Mandatory-设置目标几何形状
    seg.setModelType(pcl::SACMODEL_PLANE);
    //分割方法：随机采样法
    seg.setMethodType(pcl::SAC_RANSAC);
    //设置误差容忍范围，也就是我说过的阈值
    seg.setDistanceThreshold(0.01);

    vector<Plane> detectedLocalPlanes;

    pcl::PointCloud<pcl::PointXYZI>::Ptr curPlanePoints(new pcl::PointCloud<pcl::PointXYZI>());

    for (size_t i = 0 ; i < clusters.size() ; ++i)
    {
        Plane plane;
        pcl::PointCloud<pcl::PointXYZI>::Ptr pointsOnPlane(new pcl::PointCloud<pcl::PointXYZI>());
        pointsOnPlane->points.clear();

        for (size_t j = 0 ; j < clusters[i].indices.size() ; ++j)
        {
            pcl::PointXYZI tmpPoint;
            tmpPoint.x = laserCloudNgAftEe->points[clusters[i].indices[j]].x;
            tmpPoint.y = laserCloudNgAftEe->points[clusters[i].indices[j]].y;
            tmpPoint.z = laserCloudNgAftEe->points[clusters[i].indices[j]].z;
            tmpPoint.intensity = i;
            pointsOnPlane->points.push_back(tmpPoint);
            laserCloudNgAftReg->points.push_back(tmpPoint);
        }
        //输入点云
        seg.setInputCloud(pointsOnPlane);
        //分割点云，获得平面和法向量
        seg.segment(*inliers, *coefficients);

        double planen[4];
        planen[0] = coefficients->values[0];
        planen[1] = coefficients->values[1];
        planen[2] = coefficients->values[2];
        planen[3] = -1.0 * coefficients->values[3];

        plane.v3normal = Vector3f(planen[0], planen[1], planen[2]);
        plane.d = planen[3];
        pcl::copyPointCloud(*pointsOnPlane, *plane.planePointCloudPtr);

        pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudInlier(new pcl::PointCloud<pcl::PointXYZI>());

        for (std::vector<int>::const_iterator pit = inliers->indices.begin (); pit != inliers->indices.end (); ++pit)
        {
            laserCloudInlier->points.push_back (pointsOnPlane->points[*pit]);
        }

        pcl::copyPointCloud(*laserCloudInlier, *plane.InplanePointCloudOriPtr);


        plane.calcConvexHull(plane.planePointCloudPtr,planen);

        plane.computeMassCenterAndArea();
        plane.areaVoxels= plane.planePointCloudPtr->size() * 0.0025;

        plane.keyFrameId = keyFrameCount;

        *curPlanePoints += *plane.planePointCloudPtr;

        detectedLocalPlanes.push_back(plane);
    }

//    laserCloudOri_mp1 += *curPlanePoints;

    std::cout<<"detectedLocalPlanes.size() "<<detectedLocalPlanes.size ()<<std::endl;

    // Merge detected planes with previous ones if they are the same
    observedPlanes.clear();
    vPlanes.clear();

    for (size_t i = 0; i < detectedLocalPlanes.size (); i++)
    {
        detectedLocalPlanes[i].id = vPlanes.size();

        // Update co-visibility graph
        for(set<unsigned>::iterator it = observedPlanes.begin(); it != observedPlanes.end(); it++)
        {
            detectedLocalPlanes[i].neighborPlanes[*it] = 1;
            vPlanes[*it].neighborPlanes[detectedLocalPlanes[i].id] = 1;
        }
        observedPlanes.insert(detectedLocalPlanes[i].id);
        vPlanes.push_back(detectedLocalPlanes[i]);
    }

    for(set<unsigned>::iterator it = observedPlanes.begin(); it != observedPlanes.end(); it++)
    {
        Plane &observedPlane = vPlanes[*it];
        observedPlane.calcElongationAndPpalDir();
    }


    double time_endinSeg = pcl::getTime();

    std::cout<<"plane seg time: "<<time_endinSeg-time_behinSeg<<std::endl;

//    std::shared_ptr<Maps_keyframe<float>> mk = std::make_shared<Maps_keyframe<float>>();

    std::cout<<"vPlanes.size() ";
    std::cout<<vPlanes.size()<<std::endl;

    mk->vPlanes = vPlanes;
    mk->m_accumulated_point_cloud = *c_keyframe.orilaserCloud;

    mk->m_accumulated_ng_pc = c_keyframe.g_laserCloud;

    mk->m_keyframe_idx = keyFrameCount;
    mk->m_ending_frame_idx = c_keyframe.m_ending_frame_idx;
    mk->m_accumulate_frames = c_keyframe.framecount;
    mk->m_pose_q = c_keyframe.m_pose_q;
    mk->m_pose_t = c_keyframe.m_pose_t;
    mk->m_accumulated_line_pc = *c_keyframe.linelaserCloud;
    mk->m_accumulated_surf_pc = *c_keyframe.surflaserCloud;
//    mk->m_accumulated_g_pc = c_keyframe.g_laserCloud;

    ////////////////

    Eigen::Quaterniond                                 q_curr;
    Eigen::Vector3d                                    t_curr;

    q_curr = c_keyframe.m_pose_q;
    t_curr = c_keyframe.m_pose_t;

    keyframe_vec.push_back( mk );


    map_id_pc.insert( std::make_pair( map_id_pc.size(), keyframe_vec.back()->m_accumulated_point_cloud ) );


    int curren_frame_idx = keyframe_vec.back()->m_ending_frame_idx;

    pose3d_vec.push_back( Pose3d( q_curr, t_curr ) );
    pose3d_map.insert( std::make_pair( pose3d_map.size(), Pose3d( q_curr, t_curr ) ) );

    if ( pose3d_vec.size() >= 2 )
    {
        Constraint3d temp_csn;
        Eigen::Vector3d relative_T = pose3d_vec[ pose3d_vec.size() - 2 ].q.inverse() * ( t_curr - pose3d_vec[ pose3d_vec.size() - 2 ].p );
        Eigen::Quaterniond relative_Q = pose3d_vec[ pose3d_vec.size() - 2 ].q.inverse() * q_curr;
        temp_csn = Constraint3d( pose3d_vec.size() - 2, pose3d_vec.size() - 1,relative_Q, relative_T );
        constrain_vec.push_back( temp_csn );
    }

    std::shared_ptr<Maps_keyframe<float>>& last_keyframe = keyframe_vec.back();

    printf( "--- Current_idx = %d, lidar_frame_idx = %d ---\r\n", keyframe_vec.size(), curren_frame_idx );

    std::map<unsigned, unsigned> bestMatch; bestMatch.clear();
    int this_kf = 0;

    float wdif_height;
    float wdif_height2;
    float wdif_normal;
    float wrel_dist_centers;
    float wal;
    float wea;


    pcl::PointXYZI pointSel;
    pointSel.x = t_curr(0);
    pointSel.y = t_curr(1);
    pointSel.z = t_curr(2);

    if (all_kfs_pose.points.size() < 10) {
        all_kfs_pose.push_back(pointSel);
        return;
    }

    m_kdtree_kfs.setInputCloud( all_kfs_pose.makeShared() );
    std::vector<int> pointSearchInd;
    std::vector<float> pointSearchSqDis;
    m_kdtree_kfs.nearestKSearch(pointSel, 10, pointSearchInd, pointSearchSqDis);
    all_kfs_pose.push_back(pointSel);

    bestMatch.clear();
    this_kf = 0;

    double aveTime = 0.0;
    int graphmatchTimes = 0;
    if (last_keyframe->vPlanes.size() > 2) {

        for (size_t his = 0; his < pointSearchInd.size(); his++) {
            if ((keyframe_vec.size() - pointSearchInd[his]) <= 10) {
                continue;
            }

            if (keyframe_vec[pointSearchInd[his]]->vPlanes.size() < 3)
                continue;

            Subgraph currentSubgraph(keyframe_vec[pointSearchInd[his]]->vPlanes, 0);
            Subgraph targetSubgraph(last_keyframe->vPlanes, 0);

            int unaycount;
            double time_begComp = pcl::getTime();
            std::map<unsigned, unsigned> resultingMatch = matcher.compareSubgraphs(currentSubgraph, targetSubgraph,
                                                                                   unaycount);

            graphmatchTimes++;

            double time_endComp = pcl::getTime();

            aveTime += time_endComp - time_begComp;

            if (resultingMatch.size() > bestMatch.size()) {
                bestMatch = resultingMatch;
                this_kf = pointSearchInd[his];

                wdif_height = matcher.wdif_height;
                wdif_height2 = matcher.wdif_height2;
                wdif_normal = matcher.wdif_normal;
                wrel_dist_centers = matcher.wrel_dist_centers;
                wal = matcher.wal;
                wea = matcher.wea;
            }
        }
    }

//    std::cout<<"!!!!!!!!!"<<std::endl;
//    std::cout<<"wdif_height "<<wdif_height<<std::endl;
//    std::cout<<"wdif_height2 "<<wdif_height2<<std::endl;
//    std::cout<<"wdif_normal "<<wdif_normal<<std::endl;
//    std::cout<<"wrel_dist_centers "<<wrel_dist_centers<<std::endl;
//    std::cout<<"!!!!!!!!!"<<std::endl;


    std::cout<<"comp time "<<aveTime/graphmatchTimes<<std::endl;

    std::vector<Eigen::Vector3f> kvc;
    std::vector<Eigen::Vector3f> lvc;
    std::vector<Eigen::Vector3f> kvn;
    std::vector<Eigen::Vector3f> lvn;

    for (map<unsigned, unsigned>::iterator it = bestMatch.begin(); it != bestMatch.end(); it++)
    {
        kvc.push_back(keyframe_vec[this_kf]->vPlanes[it->first].v3center);
        lvc.push_back(last_keyframe->vPlanes[it->second].v3center);
        kvn.push_back(keyframe_vec[this_kf]->vPlanes[it->first].v3normal);
        lvn.push_back(last_keyframe->vPlanes[it->second].v3normal);
    }

    std::map<unsigned, unsigned> bestMatchVehicle; bestMatchVehicle.clear();
    int this_kf_vehicle = 0;

    if (last_keyframe->vVehicles.size() > 2) {

        for (size_t his = 0; his < pointSearchInd.size(); his++) {
            if ((keyframe_vec.size() - pointSearchInd[his]) <= 10) {
                continue;
            }

            if (keyframe_vec[pointSearchInd[his]]->vVehicles.size() < 3)
                continue;

            Subgraph currentSubgraph(keyframe_vec[pointSearchInd[his]]->vVehicles, 0);
            Subgraph targetSubgraph(last_keyframe->vVehicles, 0);

            int unaycount;
            std::map<unsigned, unsigned> resultingMatch = matcher.compareSubgraphsVehiclePlaneRef(currentSubgraph, targetSubgraph, unaycount, kvc,
                                                                                                  lvc, kvn, lvn);

            if (resultingMatch.size() > bestMatchVehicle.size()) {
                bestMatchVehicle = resultingMatch;
                this_kf_vehicle = pointSearchInd[his];
            }
        }
    }


    std::map<unsigned, unsigned> bestMatchPole; bestMatchPole.clear();
    int this_kf_pole = 0;

    if (last_keyframe->vPoles.size() > 2) {

        for (size_t his = 0; his < pointSearchInd.size(); his++) {
            if ((keyframe_vec.size() - pointSearchInd[his]) <= 10) {
                continue;
            }

            if (keyframe_vec[pointSearchInd[his]]->vPoles.size() < 3)
                continue;

            Subgraph currentSubgraph(keyframe_vec[pointSearchInd[his]]->vPoles, 0);
            Subgraph targetSubgraph(last_keyframe->vPoles, 0);

            int unaycount;
            std::map<unsigned, unsigned> resultingMatch = matcher.compareSubgraphsPolePlaneRef(currentSubgraph, targetSubgraph, unaycount, kvc,
                                                                                                  lvc, kvn, lvn);

            if (resultingMatch.size() > bestMatchPole.size()) {
                bestMatchPole = resultingMatch;
                this_kf_pole = pointSearchInd[his];
            }
        }
    }


    std::cout<<"+++++++++++++++++++++++++++++++"<<std::endl;
    std::cout<<"this_kf "<<this_kf<<std::endl;
    std::cout<<"this_kf_pole "<<this_kf_pole<<std::endl;
    std::cout<<"this_kf_vehicle "<<this_kf_vehicle<<std::endl;
    std::cout<<"=================="<<std::endl;
    std::cout<<"bestMatch.size() "<<bestMatch.size()<<std::endl;
    std::cout<<"bestMatchPole.size() "<<bestMatchPole.size()<<std::endl;
    std::cout<<"bestMatchVehicle.size() "<<bestMatchVehicle.size()<<std::endl;
    std::cout<<"-----------------------------"<<std::endl;

    if (bestMatch.size() < 3)
        return;

    int allMatchedSize = bestMatch.size();

    if (bestMatchVehicle.size() > 1 && abs(this_kf - this_kf_vehicle) < 2)
    {
        allMatchedSize += bestMatchVehicle.size();
    }


    if (bestMatchPole.size() > 1 && abs(this_kf - this_kf_pole) < 2)
    {
        allMatchedSize += bestMatchPole.size();
    }


//    std::cout<<"comp time "<<aveTime/(pointSearchInd.size())<<std::endl;

    std::cout<<"allMatchedSize "<<allMatchedSize<<std::endl;

//    if (allMatchedSize < 6) {
//
//        if (allMatchedSize < 3)
//            return;
//
//        float minarea = 0.0;
//        for (map<unsigned, unsigned>::iterator it = bestMatch.begin(); it != bestMatch.end(); it++)
//            minarea += last_keyframe->vPlanes[it->second].areaHull;
//
//        float minarea2 = 0.0;
//        for (map<unsigned, unsigned>::iterator it = bestMatch.begin(); it != bestMatch.end(); it++)
//            minarea2 += keyframe_vec[this_kf]->vPlanes[it->first].areaHull;
//
//        std::cout<<"minarea/minarea2 "<<minarea/minarea2<<std::endl;
//
//
//
//        float areaG = 0.0;
//        for (size_t his = 0; his < last_keyframe->vPlanes.size(); ++his) {
//            areaG += last_keyframe->vPlanes[his].areaHull;
//        }
//
//        float radioArea = minarea / areaG;
//
//        std::cout << "radioArea " << radioArea << std::endl;
//
////        if (radioArea < 0.30)
////            return;
//
//        std::vector<int> firstINdexes;
//        for (map<unsigned, unsigned>::iterator it = bestMatch.begin(); it != bestMatch.end(); it++) {
//            firstINdexes.push_back(it->first);
//        }
//
//        pcl::PointCloud<pcl::PointXYZI> normPOints;
//        for (int it = 0; it < firstINdexes.size(); ++it) {
//            pcl::PointXYZI phe;
//            phe.x = keyframe_vec[this_kf]->vPlanes[firstINdexes[it]].v3normal(0);
//            phe.y = keyframe_vec[this_kf]->vPlanes[firstINdexes[it]].v3normal(1);
//            phe.z = keyframe_vec[this_kf]->vPlanes[firstINdexes[it]].v3normal(2);
//            normPOints.points.push_back(phe);
//        }
//
//        //聚类
//        pcl::search::KdTree<pcl::PointXYZI>::Ptr tree1(new pcl::search::KdTree<pcl::PointXYZI>);
//        tree1->setInputCloud(normPOints.makeShared());
//        std::vector<pcl::PointIndices> cluster_indices2;
//        pcl::EuclideanClusterExtraction<pcl::PointXYZI> ec1;
//        ec1.setClusterTolerance(0.3);
//        ec1.setMinClusterSize(1);
//        ec1.setMaxClusterSize(100);
//        ec1.setSearchMethod(tree1);
//        ec1.setInputCloud(normPOints.makeShared());
//        ec1.extract(cluster_indices2);
//
//////        if (cluster_indices2.size() == 2 && cluster_indices2[0].indices.size() == 3 && cluster_indices2[1].indices.size() == 1)
//        if (cluster_indices2.size() < 2)
//            return;
//
//        if (minarea/minarea2 < 0.85 || minarea/minarea2 > 1.18)
//            return;
//    }

    if (allMatchedSize < 6) {
        return;
    }


    double icp_score = 0.0;

    double icp_score_t = 0.0;

    //////////////////////////////
    std::vector<pcl::PointCloud<pcl::PointXYZI>> v_selectedPc;
    std::vector<Eigen::Vector3d> v_selectedNorm;
    std::vector<Eigen::Vector3d> v_selectedNormMatch;

    std::vector<double> v_d;

    for (map<unsigned, unsigned>::iterator it = bestMatch.begin(); it != bestMatch.end(); it++)
    {
        pcl::PointCloud<pcl::PointXYZI> oriPlanePs;
        oriPlanePs = *keyframe_vec[this_kf]->vPlanes[it->first].InplanePointCloudOriPtr;

        pcl::VoxelGrid< pcl::PointXYZI > filterTmp;
        filterTmp.setInputCloud( oriPlanePs.makeShared() );
        filterTmp.setLeafSize( 1.0, 1.0, 1.0 );
        filterTmp.filter( oriPlanePs );
        v_selectedPc.push_back(oriPlanePs);
        Eigen::Vector3d selct_n = last_keyframe->vPlanes[it->second].v3normal.cast<double>();
        v_selectedNorm.push_back(selct_n);

        Eigen::Vector3d selct_n1 = keyframe_vec[this_kf]->vPlanes[it->first].v3normal.cast<double>();
        v_selectedNormMatch.push_back(selct_n1);

        v_d.push_back(last_keyframe->vPlanes[it->second].d);
    }

    std::vector<Eigen::Vector3d> v_centriods;
    std::vector<Eigen::Vector3d> v_centriodsMatch;

    if (bestMatchVehicle.size() > 1 && abs(this_kf - this_kf_vehicle) < 2)
    {
        for (map<unsigned, unsigned>::iterator it = bestMatchVehicle.begin(); it != bestMatchVehicle.end(); it++)
        {
            Eigen::Vector3d selct_n = last_keyframe->vVehicles[it->second].v3center.cast<double>();
            v_centriods.push_back(selct_n);

            Eigen::Vector3d selct_n1 = keyframe_vec[this_kf_vehicle]->vVehicles[it->first].v3center.cast<double>();
            v_centriodsMatch.push_back(selct_n1);
        }
    }


    if (bestMatchPole.size() > 1 && abs(this_kf - this_kf_pole) < 2)
    {
        for (map<unsigned, unsigned>::iterator it = bestMatchPole.begin(); it != bestMatchPole.end(); it++)
        {
            Eigen::Vector3d selct_n = last_keyframe->vPoles[it->second].v3center.cast<double>();
            v_centriods.push_back(selct_n);

            Eigen::Vector3d selct_n1 = keyframe_vec[this_kf_pole]->vPoles[it->first].v3center.cast<double>();
            v_centriodsMatch.push_back(selct_n1);

        }
    }

    //////////////////////////////

    Eigen::Quaterniond quaternion;
    Eigen::Vector3d    trans(0,0,0);


    double init_cost, final_cost;

    double time_begMapAlign = pcl::getTime();

    if (v_selectedNorm.size() > 2)
        scene_align.find_init_tranfrom_of_two_mappings(v_selectedNorm, v_selectedNormMatch, quaternion);

    if (v_centriods.size() > 2)
        scene_align.find_init_tranfrom_of_two_mappings2(v_centriods, v_centriodsMatch, trans);

    /************/
//    std::cout<<"====test===="<<std::endl;
//    std::cout<<keyframe_vec[this_kf]->m_accumulated_point_cloud.size()<<"  "<<keyframe_vec[this_kf]->m_accumulated_ng_pc.size()<<std::endl;
//    std::cout<<float(keyframe_vec[this_kf]->m_accumulated_point_cloud.size())/float(keyframe_vec[this_kf]->m_accumulated_ng_pc.size())<<std::endl;
//
//    scene_align.find_tranfrom_of_two_mappings(keyframe_vec[this_kf]->m_accumulated_point_cloud, last_keyframe->m_accumulated_point_cloud, icp_score_t);
//    scene_align.find_tranfrom_of_two_mappings(keyframe_vec[this_kf]->m_accumulated_ng_pc, last_keyframe->m_accumulated_ng_pc, icp_score_t);
//    std::cout<<"====test===="<<std::endl;
    /************/

    scene_align.find_tranfrom_of_two_mappings(keyframe_vec[this_kf]->m_accumulated_surf_pc, keyframe_vec[this_kf]->m_accumulated_line_pc,
                                              last_keyframe->m_accumulated_surf_pc, last_keyframe->m_accumulated_line_pc,
                                              icp_score, v_selectedPc, v_selectedNorm, v_d, quaternion, trans, init_cost, final_cost);

    double time_endMapAlign = pcl::getTime();

    std::cout<<"align time "<<time_endMapAlign-time_begMapAlign<<std::endl;


    std::cout<<icp_score<<"  "<<init_cost<<"  "<<final_cost<<std::endl;

    std::cout<<"=============**************"<<std::endl;

    printf("ICP inlier threshold = %lf\r\n", icp_score);

    auto Q_a = pose3d_vec[this_kf].q;
    auto Q_b = pose3d_vec[pose3d_vec.size() - 1].q;
    auto T_a = pose3d_vec[this_kf].p;
    auto T_b = pose3d_vec[pose3d_vec.size() - 1].p;
    auto ICP_q = scene_align.m_q_w_curr;
    auto ICP_t = scene_align.m_t_w_curr;

    ICP_t = ( ICP_q.inverse() * ( -ICP_t ) );
    ICP_q = ICP_q.inverse();

    std::cout << "ICP_q = " << ICP_q.coeffs().transpose() << std::endl;
    std::cout << "ICP_t = " << ICP_t.transpose() << std::endl;

    if ( icp_score < 0.31 )

//    if ( icp_score < 0.31 && final_cost > 50.0)
    {

        v_icp.push_back(icp_score_t);

        double t_s = 0.0;
        for (auto item : v_icp)
        {
            t_s += item;
        }

        std::cout<<"t_s ave: "<<t_s/v_icp.size()<<std::endl;

//        for (map<unsigned, unsigned>::iterator it = bestMatch.begin(); it != bestMatch.end(); it++)
//        {
//            laserCloudOri_m2 += *keyframe_vec[this_kf]->vPlanes[it->first].planePointCloudPtr;
//            laserCloudOri_m2 += *last_keyframe->vPlanes[it->second].planePointCloudPtr;
//        }

        for (map<unsigned, unsigned>::iterator it = bestMatchVehicle.begin(); it != bestMatchVehicle.end(); it++)
        {
            laserCloudOri_m2 += *keyframe_vec[this_kf_vehicle]->vVehicles[it->first].VehiclePointCloudPtr;
            laserCloudOri_m2 += *last_keyframe->vVehicles[it->second].VehiclePointCloudPtr;
        }

        std::cout<<"bestMatchVehicle.size() "<<bestMatchVehicle.size()<<std::endl;
        std::cout<<"last_keyframe->vVehicles.size() "<<last_keyframe->vVehicles.size()<<std::endl;
        std::cout<<"keyframe_vec[this_kf_vehicle]->vVehicles.size() "<<keyframe_vec[this_kf_vehicle]->vVehicles.size()<<std::endl;


        printf("I believe this is true loop.\r\n");

        VectorOfConstraints constrain_vec_temp;
        constrain_vec_temp = constrain_vec;

        /*闭环约束*/
        Constraint3d pose_constrain;
        auto q_res = Q_b.inverse() * ICP_q.inverse() * Q_a;
        //q_res = q_res.inverse();
        auto t_res = Q_b.inverse() * ( ICP_q.inverse() * ( T_a - ICP_t ) - T_b );
        //t_res = q_res.inverse()*(-t_res);
        //q_res = q_res.inverse();

//                cout << "=== Add_constrain_of_loop ====" << endl;
//                cout << Q_a.coeffs().transpose() << endl;
//                cout << Q_b.coeffs().transpose() << endl;
//                cout << ICP_q.coeffs().transpose() << endl;
//                cout << T_a.transpose() << endl;
//                cout << T_b.transpose() << endl;
//                cout << ICP_t.transpose() << endl;
//                cout << "Result: " << endl;
//                cout << q_res.coeffs().transpose() << endl;
//                cout << t_res.transpose() << endl;

        //t_res.setZero();
        pose_constrain.id_begin = pose3d_vec.size() - 1;
        pose_constrain.id_end = this_kf;
        pose_constrain.t_be.p = t_res;
        pose_constrain.t_be.q = q_res;
        /*闭环约束*/


        loop_closure_matchedid[pose3d_vec.size() - 1] = this_kf;

        constrain_vec_temp.push_back( pose_constrain );


        pose3d_map_ori = pose3d_map;
        auto temp_pose_3d_map = pose3d_map;

        pose_graph_optimization(temp_pose_3d_map, constrain_vec_temp);

        optimized_pose3d_map = temp_pose_3d_map;

        ////////

//        constrain_vec = constrain_vec_temp;
//        pose3d_map = temp_pose_3d_map;

        ////////

        OutputPoses( std::string("/home/jjwen/result/poses_ori.txt" ), pose3d_map_ori );
        OutputPoses( std::string( "/home/jjwen/result/poses_opm.txt" ), temp_pose_3d_map );


//        refined_pt = map_rfn.refine_pointcloud( map_id_pc, pose3d_map_ori, temp_pose_3d_map, ( int ) map_id_pc.size() - 1);

        refined_pt.points.clear();
        refined_pt_bef.points.clear();

        for ( int pc_idx = ( int ) map_id_pc.size() - 1; pc_idx >= 0; pc_idx -= 2 )
        {
//            std::cout << "*** Refine pointcloud, curren idx = " << pc_idx << " ***" << endl;
            refined_pt += map_rfn.refine_pointcloud( map_id_pc, pose3d_map_ori, temp_pose_3d_map, pc_idx);

            refined_pt_bef += map_rfn.refine_pointcloud( map_id_pc, pose3d_map_ori, pose3d_map_ori, pc_idx);


        }
        //map_rfn.refine_mapping( path_name, 0 );
        if ( 0 )
        {
            map_rfn.refine_mapping( map_id_pc, pose3d_map_ori, temp_pose_3d_map);
//            map_rfn.m_pts_aft_refind
        }

    }

}





void LiPMatch::run()
{
    size_t numPrevKFs = 0;

    while(!LiPMatch_stop)  // Stop loop if LiPMatch
    {
        if( numPrevKFs == frameQueue.size() )
        {
          sleep(10);
        }
        else
        {
            detectPlanesCloud( frameQueue[numPrevKFs], numPrevKFs);
            ++numPrevKFs;
        }
    }
    LiPMatch_finished = true;
}



bool LiPMatch::stop_LiPMatch()
{
    LiPMatch_stop = true;
    while(!LiPMatch_finished)
        sleep(1);

    cout << "Waiting for LiPMatch thread to die.." << endl;

    joinThread(LiPMatch_hd);
    LiPMatch_hd.clear();

    return true;
}

LiPMatch::~LiPMatch()
{
    cout << "\n\n\nLiPMatch destructor called -> Save color information to file\n";

    stop_LiPMatch();

    cout << " .. LiPMatch has died." << endl;
}


