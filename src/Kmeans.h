/*
 * Kmeans.h
 *
 *  Created on: 04-Mar-2009
 *      Author: sbutler
 */

#pragma once

//#include <Eigen/Core>
//#include <Eigen/Array>
//#include <vector>

#include <igl/principal_curvature.h>


class Kmeans {
public:
	Kmeans();
	virtual ~Kmeans();

	/**
	 * kmeans clustering
	 * based on kmeans matlab script by Ian T Nabney (1996-2001)
	 *
	 * @param data 	affinity matrix
	 * @oaram ncentres	number of clusters
	 * @return		ordered set of data row indices (the path id) for each cluster
	 * 				(order is based on distance to cluster centre)
	 */
	static std::vector<std::vector<int> > cluster(Eigen::MatrixXd& data, int ncentres);
};

