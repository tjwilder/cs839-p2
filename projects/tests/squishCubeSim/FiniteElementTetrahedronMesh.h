#pragma once

#include "AnimatedTetrahedronMesh.h"

#include "CGVectorWrapper.h"
#include "CGSystemWrapper.h"
#include "ConjugateGradient.h"

#include <Eigen/Dense>

#define USE_LINEAR_ELASTICITY
//#define USE_ST_VENANT_KIRCHHOFF

template<class T>
struct FiniteElementTetrahedronMesh : public AnimatedTetrahedronMesh<T>
{
    using Base = AnimatedTetrahedronMesh<T>;
    using Base::m_meshElements;
    using Base::m_particleX;
    using Vector3 = typename Base::Vector3;
    using Matrix33 = Eigen::Matrix< T , 3 , 3>;

    int m_nFrames;
    int m_subSteps;
    T m_frameDt;
    T m_stepDt;
    T m_stepEndTime;

    const T m_density;
    const T m_mu;
    const T m_lambda;
    const T m_rayleighCoefficient;
    const T m_singularValueThreshold;
    
    std::vector<T> m_particleMass;
    std::vector<Matrix33> m_DmInverse;
    std::vector<T> m_restVolume;
	mutable std::vector<Matrix33> m_F;
    
	FiniteElementTetrahedronMesh(const T density, const T mu, const T lambda, const T rayleighCoefficient)
        :m_density(density), m_mu(mu), m_lambda(lambda), m_rayleighCoefficient(rayleighCoefficient), m_singularValueThreshold(.2f)
    {}

    void initializeUndeformedConfiguration()
    {
        // Initialize rest shape and particle mass (based on constant density)
        m_particleMass.resize(m_particleX.size(), T()); // Initialize all particle masses to zero
        for(const auto& element: m_meshElements)
        {
            Matrix33 Dm;
            for(int j = 0; j < 3; j++)
                Dm.col(j) = m_particleX[element[j+1]] - m_particleX[element[0]];

            T restVolume = .5 * Dm.determinant();
            if(restVolume < 0)
                throw std::logic_error("Inverted element");
            m_DmInverse.emplace_back(Dm.inverse());
            m_restVolume.push_back(restVolume);
            T elementMass = m_density * restVolume;
            for(const int v: element)
                m_particleMass[v] += (1./4.) * elementMass;
        }
    }
    
    void addElasticForce(std::vector<Vector3>& f) const
    {
		m_F.clear();
        for(int e = 0; e < m_meshElements.size(); e++)
        {
            const auto& element = m_meshElements[e];

            // Compute deformation gradient
            Matrix33 Ds = Vector3::Zero().asDiagonal();
			for (int j = 0; j < 3; j++) {
				Ds.col(j) = m_particleX[element[j + 1]] - m_particleX[element[0]];
			}
            Matrix33 F = Ds * m_DmInverse[e];

            // Compute SVD
			Eigen::JacobiSVD<Matrix33> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Matrix33 U = svd.matrixU();
            Matrix33 V = svd.matrixV();
            Vector3 vSigma = svd.singularValues();
            if ( U.determinant() < 0. ) {
                if ( V.determinant() < 0. ) {
                    // Both determinants negative, just negate 3rd column on both
                    U.col(2) *= -1.f;
                    V.col(2) *= -1.f;
                }
                else {
                    // Only U has negative determinant, negate 3rd column and second singular value
                    U.col(2) *= -1.f;
                    vSigma[2] = -vSigma[2];
                }
            }
			else {
				if (V.determinant() < 0.) {
					// Only V has negative determinant, negate 3rd column and second singular value
					V.col(2) *= -1.f;
					vSigma[2] = -vSigma[2];
				}
			}
            if ( (F-U*vSigma.asDiagonal()*V.transpose()).norm() > 1e-5 )
                throw std::logic_error("SVD error");

            // Apply thresholding of singular values, and re-constitute F
            for (int v = 0; v < 3; v++)
                vSigma[v] = std::max<T>(m_singularValueThreshold, vSigma[v]);
            Matrix33 Sigma = vSigma.asDiagonal();
            F = U * Sigma * V.transpose();

			m_F.push_back(F);
            
#ifdef USE_LINEAR_ELASTICITY
            Matrix33 strain = .5 * (F + F.transpose()) - Matrix33::Identity();
            Matrix33 P = 2. * m_mu * strain + m_lambda * strain.trace() * Matrix33::Identity();
#endif

#ifdef USE_ST_VENANT_KIRCHHOFF
            Matrix33 E = .5 * ( F.transpose() * F - Matrix33::Identity());
            Matrix33 P = F * (2. * m_mu * E + m_lambda * E.trace() * Matrix33::Identity());
#endif

            Matrix33 H = -m_restVolume[e] * P * m_DmInverse[e].transpose();
            
            for(int j = 0; j < 3; j++){
                f[element[j+1]] += H.col(j);
                f[element[0]] -= H.col(j);
            }
        }
    }

    void addProductWithStiffnessMatrix(std::vector<Vector3>& dx, std::vector<Vector3>& df, const T scale) const
    {
        for(int e = 0; e < m_meshElements.size(); e++)
        {
            const auto& element = m_meshElements[e];

            // Compute deformation gradient
            Matrix33 Ds;
            for(int j = 0; j < 3; j++)
                Ds.col(j) = m_particleX[element[j+1]] - m_particleX[element[0]];

			Matrix33 F = m_F[e];
            
            // Compute differential(s)
            Matrix33 dDs;
            for(int j = 0; j < 3; j++)
                dDs.col(j) = dx[element[j+1]] - dx[element[0]];
            Matrix33 dF = dDs * m_DmInverse[e];

#ifdef USE_LINEAR_ELASTICITY
            Matrix33 dstrain = .5 * (dF + dF.transpose());
            Matrix33 dP = scale * (2. * m_mu * dstrain + m_lambda * dstrain.trace() * Matrix33::Identity());
#endif

#ifdef USE_ST_VENANT_KIRCHHOFF
            Matrix33 E = .5 * ( F.transpose() * F - Matrix33::Identity());
            Matrix33 dE = .5 * ( dF.transpose() * F + F.transpose() * dF);
            Matrix33 dP = dF * (2. * m_mu *  E + m_lambda *  E.trace() * Matrix33::Identity()) +
                           F * (2. * m_mu * dE + m_lambda * dE.trace() * Matrix33::Identity());

#endif
            
            Matrix33 dH = m_restVolume[e] * dP * m_DmInverse[e].transpose();
            
            for(int j = 0; j < 3; j++){
                df[element[j+1]] += dH.col(j);
                df[element[0]] -= dH.col(j);
            }
        }
    }

    void simulateSubstep()
    {
        using FEMType = FiniteElementTetrahedronMesh<T>;        

        const int nParticles = m_particleX.size();

        // Construct initial guess for next-timestep
        //   Positions -> Same as last timestep
        
        // Overwrite boundary conditions with desired values

        setBoundaryConditions();
        
        // Solve for everything else using Conjugate Gradients

        std::vector<Vector3> dx(nParticles, Vector3::Zero());
        std::vector<Vector3> rhs(nParticles, Vector3::Zero());
        std::vector<Vector3> q(nParticles, Vector3::Zero());
        std::vector<Vector3> s(nParticles, Vector3::Zero());
        std::vector<Vector3> r(nParticles, Vector3::Zero());
        CGVectorWrapper<Vector3> dxWrapper(dx);
        CGVectorWrapper<Vector3> rhsWrapper(rhs);
        CGVectorWrapper<Vector3> qWrapper(q);
        CGVectorWrapper<Vector3> sWrapper(s);
        CGVectorWrapper<Vector3> rWrapper(r);
        CGSystemWrapper<Vector3, FEMType> systemWrapper(*this);
        
        addElasticForce(rhs);
        clearConstrainedParticles(rhs);

        ConjugateGradient<T>::Solve(systemWrapper,
            dxWrapper, rhsWrapper, qWrapper, sWrapper, rWrapper,
            1e-4, 50);

        // Apply corrections to positions and velocities

        for(int p = 0; p < nParticles; p++)
            m_particleX[p] += dx[p];
    }

    void simulateFrame(const int frame)
    {
        m_stepDt = m_frameDt / (T) m_subSteps;

		std::cout << "Simulating frame " << frame << std::endl;
        for(int step = 1; step <= m_subSteps; step++){
            m_stepEndTime = m_frameDt * (T) (frame-1) + m_stepDt * (T) step;
            simulateSubstep();
        }
    }

    virtual void clearConstrainedParticles(std::vector<Vector3>& x) {}
    virtual void setBoundaryConditions() {}
};