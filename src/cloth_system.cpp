#include <cloth_system.hpp>
#include <shader_c.hpp>
#include <shader_s.hpp>
#include <gl_bindings.h>
#include <chrono>
#include <iostream>
#include<Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include<Eigen/SparseCholesky>
#include<Eigen/IterativeLinearSolvers>


#if PARDISO_SOLVE == 1
#include <pardiso.h>
#endif

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include "stb_image.h"


#define EIGEN_DONT_PARALLELIZE



static double signedAngle(
    const Eigen::Ref<const Eigen::Vector3d>& u, 
    const Eigen::Ref<const Eigen::Vector3d>& v, 
    const Eigen::Ref<const Eigen::Vector3d>& n){
    Eigen::Vector3d w = u.cross(v);
    double angle = atan2(w.norm(), u.dot(v));
    if (n.dot(w) < 0){
        angle = -angle;
    }
    return angle;
}


static Eigen::Matrix3d mmt(
    const Eigen::Ref<const Eigen::Matrix3d>& matrix){
    return matrix + matrix.transpose();
}

/*
"""# Hinge angle, its gradient, and Hessian"""

#          x2
#          /\
#         /  \
#      e1/    \e3
#       /  t0  \
#      /        \
#     /    e0    \
#   x0------------x1
#     \          /
#      \   t1   /
#       \      /
#      e2\    /e4
#         \  /
#          \/
#          x3
#
#  Edge orientation: e0,e1,e2 point away from x0
#                       e3,e4 point away from x1
*/

static double getTheta(
        const Eigen::Ref< const Eigen::Vector3d > & node0, 
        const Eigen::Ref< const Eigen::Vector3d > & node1, 
        const Eigen::Ref< const Eigen::Vector3d > & node2, 
        const Eigen::Ref< const Eigen::Vector3d > & node3){
    Eigen::Vector3d m_e0 = node1 - node0;
    Eigen::Vector3d m_e1 = node2 - node0;
    Eigen::Vector3d m_e2 = node3 - node0;

    Eigen::Vector3d n0 = m_e0.cross(m_e1);
    Eigen::Vector3d n1 = m_e2.cross(m_e0);

    return signedAngle(n0, n1, m_e0);
}

void ClothSystem::InitializeSystemData(void* params) {
    ClothSystemParameters* clothParams = static_cast<ClothSystemParameters*>(params);
    clothParams->particleCount = (clothParams->clothSideLength + 1) * (clothParams->clothSideLength + 1);
    particleSystemDataBlock.particleCount = clothParams->particleCount;
    particleCount = clothParams->particleCount;
    particleSystemDataBlock.timeStep = clothParams->timeStep;
    dofCount = 3 * particleCount;
    vecValCount = particleCount * 4;
    
    indiceCount = 6 * clothParams->clothSideLength * clothParams->clothSideLength;

    double Y = clothParams->youngModulus;
    double h = clothParams->thickness;
    
    particleSystemDataBlock.bendingStiffness = 2.0 / sqrt(3.0) * Y * pow(h, 3.0) / 12.0;
    particleSystemDataBlock.timeStep = clothParams->timeStep;
    particleSystemDataBlock.tolerance = particleSystemDataBlock.bendingStiffness / (0.01) * 1.0e-3;
    particleSystemDataBlock.gravityAccel = clothParams->gravityAccel;


    bendingStiffness = particleSystemDataBlock.bendingStiffness;    
    dt = particleSystemDataBlock.timeStep;
    tol = particleSystemDataBlock.tolerance;
    gravityAccel = particleSystemDataBlock.gravityAccel;
    totalTime = 0.0;
    type = "base";
};

void ClothSystem::InitializeBufferData(void* params) {
    ClothSystemParameters* clothParams = static_cast<ClothSystemParameters*>(params);
    clothParams->particleCount = (clothParams->clothSideLength+ 1) * (clothParams->clothSideLength + 1);
    clothPositionVec.resize(vecValCount);
    particleDataVec.resize(clothParams->particleCount);
    clothTexcoords.resize(clothParams->particleCount);

    lastdofPositions.resize(dofCount);
    dofPositions.resize(dofCount);
    dofVelocities.resize(dofCount);

    glm::dvec4 meshStartPosition = glm::dvec4(-4.36647, 0.2, -1.0, 1.0);
    double yIncrement = 0.5 * sqrt(3.0) * clothParams->cellSideLength;

    uint32_t clothSideParticleCount = clothParams->clothSideLength + 1;
    double deltaMass = clothParams->totalMass / clothParams->particleCount;
    massVector.resize(dofCount, deltaMass);

    uint32_t clothSideLength = clothParams->clothSideLength;
    // All 4 corners
     fixedNodes.push_back(0);
    //  fixedNodes.push_back(clothSideLength);
     fixedNodes.push_back(clothSideParticleCount * clothSideLength);
    //  fixedNodes.push_back(clothSideParticleCount * clothSideLength + clothSideLength);

    // Round table
     //fixedNodes.push_back(clothSideLength / 2);
     //fixedNodes.push_back(clothSideParticleCount * clothSideLength / 2);
     //fixedNodes.push_back(clothSideParticleCount * clothSideLength / 2 + clothSideLength);
     //fixedNodes.push_back(clothSideParticleCount * clothSideLength + clothSideLength / 2);

    // TableCloth 
    //fixedNodes.push_back(clothSideParticleCount * (clothSideLength / 3) + (clothSideLength / 3));
    //fixedNodes.push_back(clothSideParticleCount * (clothSideLength / 3) + (2 * clothSideLength / 3));
    //fixedNodes.push_back(clothSideParticleCount * (2 * clothSideLength / 3) + (clothSideLength / 3));
    //fixedNodes.push_back(clothSideParticleCount * (2 * clothSideLength / 3) + (2 * clothSideLength / 3));

    // Handkerchief
    //fixedNodes.push_back(clothSideParticleCount * (clothSideLength / 2) + (clothSideLength / 2));
    //fixedNodes.push_back(clothSideParticleCount * (clothSideLength / 2) + (clothSideLength / 2) + 1);


    std::sort(fixedNodes.begin(), fixedNodes.end());

    uint32_t currFixedIdx = 0;
    for (uint32_t node = 0; node < particleCount; ++node) {
        if (currFixedIdx < fixedNodes.size() && node == fixedNodes[currFixedIdx]){
            ++currFixedIdx;
            fixedIndices.push_back(3 * node);
            fixedIndices.push_back(3 * node + 1);
            fixedIndices.push_back(3 * node + 2);

            oldToNewIndiceMapping.push_back(UINT32_MAX);
            oldToNewIndiceMapping.push_back(UINT32_MAX);
            oldToNewIndiceMapping.push_back(UINT32_MAX);

        
        }
        else {
            freeIndices.push_back(3 * node);
            freeIndices.push_back(3 * node + 1);
            freeIndices.push_back(3 * node + 2);

            oldToNewIndiceMapping.push_back(3 * node     - currFixedIdx * 3);
            oldToNewIndiceMapping.push_back(3 * node + 1 - currFixedIdx * 3);
            oldToNewIndiceMapping.push_back(3 * node + 2 - currFixedIdx * 3);
        }
    }
    freeDOFCount = freeIndices.size();
    fixedNodesCount = fixedNodes.size();

    massMatrix = Eigen::SparseMatrix<double,
    #if PARDISO_SOLVE == 1
        Eigen::RowMajor
    #else
        Eigen::ColMajor
    #endif
     >(freeDOFCount, freeDOFCount);
    SparseEntries massEntries;
    massEntries.reserve(freeDOFCount);
    for (uint32_t i = 0; i < freeDOFCount; ++i) {
        massEntries.push_back(Eigen::Triplet<double>(i, i, deltaMass / (dt * dt)));
    }
    massMatrix.setFromTriplets(massEntries.begin(), massEntries.end());

    // Texcoord calculations
    float maxClothWidth = (clothParams->clothSideLength + 0.5f) * clothParams->cellSideLength;
    float maxClothHeight = clothParams->clothSideLength * yIncrement;



    for (uint32_t rowIdx = 0; rowIdx < clothSideParticleCount; ++rowIdx) {
        glm::dvec4 rowStartPosition =
            meshStartPosition;

        rowStartPosition.z -= yIncrement * rowIdx;
        if (rowIdx % 2 == 1) {
            rowStartPosition.x += 0.5 * clothParams->cellSideLength;
        }

        glm::dvec4 currentPosition = rowStartPosition;
        for (uint32_t colIdx = 0; colIdx < clothSideParticleCount; ++colIdx) {
            uint32_t particleIdx = rowIdx * clothSideParticleCount + colIdx;
            for (uint32_t i = 0; i < 4; ++i) {
                clothPositionVec[particleIdx * 4 + i] = currentPosition[i];
            }
            for (uint32_t i = 0; i < 3; ++i) {
                lastdofPositions[particleIdx * 3 + i] = currentPosition[i];
                dofPositions[particleIdx * 3 + i] = currentPosition[i];
                dofVelocities[particleIdx * 3 + i] = 0.0;
            }
            clothTexcoords[particleIdx] = glm::vec2(currentPosition.x / maxClothWidth, -currentPosition.z / maxClothHeight);
            particleDataVec[particleIdx].mass = deltaMass;
            currentPosition.x += clothParams->cellSideLength;
        }
    }

    // Fill in the EBO, edges, and hinges
    uint32_t c = clothParams->clothSideLength;
    edgeCount = c * (3 * c + 2);
    hingeCount = c * (3 * c - 2);
    edgeVec.resize(edgeCount);
    hingeVec.resize(hingeCount);
    uint32_t edgeIdx = 0;
    uint32_t hingeIdx = 0;
    
    indicesVec.resize(6 * clothParams->clothSideLength * clothParams->clothSideLength);
    for (uint32_t rowIdx = 0; rowIdx < clothParams->clothSideLength; ++rowIdx) {
        if (rowIdx % 2 == 1) {
            edgeVec[edgeIdx] = glm::ivec2(rowIdx * clothSideParticleCount, (rowIdx + 1) * clothSideParticleCount);
            edgeIdx += 1;
        }   
        for (uint32_t colIdx = 0; colIdx < clothParams->clothSideLength; ++colIdx) {
            // ADC / ABD
            //    G---------C---------D---------E   
            //     \   ____/ \   ____/ \   ____/ \  
            //      \ /       \ /       \ /       \ 
            //       H---------A---------B---------F
            //
            //    ABC / BDC
            //       G---------C---------D---------E
            //      / \____   / \____   / \____   / 
            //     /       \ /       \ /       \ /  
            //    H---------A---------B---------F   

            //                          Edges = Horizontals + Verticals
            //    0---1---2             Edges = c * (c + 1) + (2 * c + 1) * c 
            //    |\  |\  |             c * [(c + 1) + (2 * c + 1)] = c * (3 * c + 2) 
            //    |  \|  \|             2 * 3 + 5 * 2 = 6 + 10 = 16
            //    3---4---5             Hinges = Diagonals + Verticals + Horizontals
            //    |  /|  /|               
            //    |/  |/  |             Hinges = c * c + (c - 1) * c + (c - 1) * c
            //    6---7---8             Hinges = c * ( 3 * c - 2 )
            //         
            

            uint32_t cIdx = rowIdx * clothSideParticleCount + colIdx;
            uint32_t dIdx = cIdx + 1;
            uint32_t aIdx = (rowIdx + 1) * clothSideParticleCount + colIdx;
            uint32_t bIdx = aIdx + 1;

            uint32_t indiceIdx = 6 * (rowIdx * clothParams->clothSideLength + colIdx);

            uint32_t topHingeIdx = cIdx - clothSideParticleCount;
            uint32_t bottomHingeIdx = aIdx;
            uint32_t rightHingeIdx = cIdx + 1;
            uint32_t leftHingeIdx = aIdx - 1;

            if (rowIdx % 2 == 0) {
                indicesVec[indiceIdx + 0] = aIdx;
                indicesVec[indiceIdx + 1] = dIdx;
                indicesVec[indiceIdx + 2] = cIdx;

                indicesVec[indiceIdx + 3] = aIdx;
                indicesVec[indiceIdx + 4] = bIdx;
                indicesVec[indiceIdx + 5] = dIdx;

                hingeVec[hingeIdx] = glm::ivec4(aIdx, dIdx, cIdx, bIdx);

                edgeVec[edgeIdx + 0] = glm::ivec2(aIdx, dIdx);
                edgeVec[edgeIdx + 1] = glm::ivec2(dIdx, cIdx);
                edgeVec[edgeIdx + 2] = glm::ivec2(cIdx, aIdx);
                
            }
            else {
                indicesVec[indiceIdx + 0] = aIdx;
                indicesVec[indiceIdx + 1] = bIdx;
                indicesVec[indiceIdx + 2] = cIdx;

                indicesVec[indiceIdx + 3] = bIdx;
                indicesVec[indiceIdx + 4] = dIdx;
                indicesVec[indiceIdx + 5] = cIdx;
                
                hingeVec[hingeIdx] = glm::ivec4(cIdx, bIdx, dIdx, aIdx);
            
                edgeVec[edgeIdx + 0] = glm::ivec2(bIdx, dIdx);
                edgeVec[edgeIdx + 1] = glm::ivec2(dIdx, cIdx);
                edgeVec[edgeIdx + 2] = glm::ivec2(cIdx, bIdx);

                topHingeIdx = dIdx - clothSideParticleCount;
                bottomHingeIdx = bIdx;
                rightHingeIdx = aIdx + 1;
                leftHingeIdx = cIdx - 1;
            }
            hingeIdx += 1;
            edgeIdx += 3;

            if (colIdx > 0){
                hingeVec[hingeIdx] = glm::ivec4(cIdx, aIdx, rightHingeIdx, leftHingeIdx);
                hingeIdx += 1;
            } 

            if (rowIdx > 0){
                hingeVec[hingeIdx] = glm::ivec4(cIdx, dIdx, topHingeIdx, bottomHingeIdx);
                hingeIdx += 1;
            }

        }

        uint32_t lastRowParticleIdx = rowIdx * clothSideParticleCount + clothParams->clothSideLength;
        if (rowIdx % 2 == 0) {
            edgeVec[edgeIdx] = glm::ivec2(lastRowParticleIdx, lastRowParticleIdx + clothSideParticleCount);
            edgeIdx += 1;
        }   


    }
    uint32_t cIdx = clothParams->clothSideLength * clothSideParticleCount;
    for (uint32_t colIdx = 0; colIdx < clothParams->clothSideLength; ++colIdx){
        edgeVec[edgeIdx] = glm::ivec2(cIdx, cIdx + 1);
        edgeIdx += 1;
        cIdx += 1;
    }

    undeformedEdgeLengthVec.resize(edgeCount);
    elasticStretchingVec.resize(edgeCount);
    
    const double Y = clothParams->youngModulus;
    const double h = clothParams->thickness;
    std::cout << std::endl;
    for (uint32_t edgeIdx = 0; edgeIdx < edgeCount; ++edgeIdx) {
        glm::ivec2 edge = edgeVec[edgeIdx];
        glm::dvec4 edgeStart (clothPositionVec[4 * edge[0]], clothPositionVec[4 * edge[0] + 1], clothPositionVec[4 * edge[0] + 2], 1.0);
        glm::dvec4 edgeEnd (clothPositionVec[4 * edge[1]], clothPositionVec[4 * edge[1] + 1], clothPositionVec[4 * edge[1] + 2], 1.0);
        glm::dvec4 edgeVec = edgeEnd - edgeStart;
        undeformedEdgeLengthVec[edgeIdx] = glm::length(edgeVec);
        elasticStretchingVec[edgeIdx] = 0.5 * sqrt(3.0) * Y * h * pow(undeformedEdgeLengthVec[edgeIdx], 2.0);       
    }

    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> bufferStrides(4, 1);
    Eigen::Map<
        Eigen::Matrix<
                double,
                Eigen::Dynamic, 1>
    > q(dofPositions.data(), dofCount, 1);
    std::cout << std::endl;
    thetaBarVec.resize(hingeCount);
    for (uint32_t hingeIdx = 0; hingeIdx < hingeCount; hingeIdx++){
        glm::ivec4 hinge = hingeVec[hingeIdx];
        Eigen::Vector3d x0 = q.segment<3>(3 * hinge[0]);  
        Eigen::Vector3d x1 = q.segment<3>(3 * hinge[1]);
        Eigen::Vector3d x2 = q.segment<3>(3 * hinge[2]);
        Eigen::Vector3d x3 = q.segment<3>(3 * hinge[3]);

        thetaBarVec[hingeIdx] = getTheta(x0, x1, x2, x3);
    }

    externalForcesVec.resize(dofCount);
    for (uint32_t particleIdx = 0; particleIdx < particleCount; ++particleIdx) {
        externalForcesVec[3 * particleIdx] = -gravityAccel * particleDataVec[particleIdx].mass;
        externalForcesVec[3 * particleIdx + 1] = - 0.5 * gravityAccel * particleDataVec[particleIdx].mass;
        externalForcesVec[3 * particleIdx + 2] = 0.0;
    }

    fmt::print("Cloth Side Length: {}\n", clothParams->clothSideLength);
    fmt::print("Total Nodes: {}\n", particleCount);
    fmt::print("Total DOF: {}\n", freeDOFCount);
    fmt::print("Total Edges: {}\n", edgeCount);
    fmt::print("Total Hinges: {}\n", hingeCount);
}

void ClothSystem::InitializeShaders(void* params) {
    ClothSystemParameters* clothParams = static_cast<ClothSystemParameters*>(params);
    std::string vertPath = "./resources/shader/" + clothParams->shaderName + "_particle.vert";
    std::string fragPath = "./resources/shader/" + clothParams->shaderName + "_particle.frag";
    std::string compPath = "./resources/shader/" + clothParams->shaderName + "_particle.comp";
    std::string edgeDebugVertPath = "./resources/shader/" + clothParams->shaderName + "_edge.vert";
    std::string edgeDebugFragPath = "./resources/shader/" + clothParams->shaderName + "_edge.frag";
    std::string hingeDebugVertPath = "./resources/shader/" + clothParams->shaderName + "_hinge.vert";
    std::string hingeDebugFragPath = "./resources/shader/" + clothParams->shaderName + "_hinge.frag";
    std::string fixedNodeVertPath = "./resources/shader/" + clothParams->shaderName + "_fixed.vert";
    std::string fixedNodeFragPath = "./resources/shader/" + clothParams->shaderName + "_fixed.frag";

    particleShader = new Shader(vertPath.c_str(), fragPath.c_str());
    particleComputeShader = new ComputeShader(compPath.c_str());

    edgeDebugShader = new Shader(edgeDebugVertPath.c_str(), edgeDebugFragPath.c_str());
    hingeDebugShader = new Shader(hingeDebugVertPath.c_str(), hingeDebugFragPath.c_str());
    fixedNodesDebugShader = new Shader(fixedNodeVertPath.c_str(), fixedNodeFragPath.c_str());
}

void ClothSystem::Clear(){
    fixedNodes.clear();
    oldToNewIndiceMapping.clear();
    freeIndices.clear();
    fixedIndices.clear();
    massVector.clear();
    massMatrix.resize(0, 0);
    dofPositions.clear();
    dofVelocities.clear();
    lastdofPositions.clear();
    externalForcesVec.clear();
    normalsVec.clear();
    clothPositionVec.clear();
    clothVelocityVec.clear();
    clothForceVec.clear();
    edgeVec.clear();
    undeformedEdgeLengthVec.clear();
    elasticStretchingVec.clear();
    thetaBarVec.clear();
    hingeVec.clear();
    indicesVec.clear();
    particleDataVec.clear();
    clothTexcoords.clear();
    

    UnbindBuffers();
    glDeleteVertexArrays(1, &fixedNodesVAO);
    
    glDeleteBuffers(1, &positionBuffer);
    glDeleteBuffers(1, &EBO);
    glDeleteBuffers(1, &edgesBuffer);
    glDeleteBuffers(1, &hingesBuffer);
    glDeleteBuffers(1, &edgeDebugDataBuffer);
    glDeleteBuffers(1, &edgeDebugEBO);
    glDeleteBuffers(1, &pickedDebugDataBuffer);
    glDeleteBuffers(1, &fixedNodesDataBuffer);
    glDeleteBuffers(1, &clothTexcoordsBuffer);
    glDeleteTextures(1, &clothTexture);
}

void ClothSystem::Reinitialize(void* params){
    Clear();
    Initialize(params);
    totalForceTime = 0.0;
    totalSolveTime = 0.0;
}

void ClothSystem::InitializeBuffers(){
    glCreateBuffers(1, &positionBuffer);
    glNamedBufferStorage(
        positionBuffer, 
        sizeof(float) * vecValCount, 
        clothPositionVec.data(),
         GL_DYNAMIC_STORAGE_BIT
    );

    glCreateBuffers(1, &EBO);
    glNamedBufferStorage(
        EBO, 
        sizeof(uint32_t) * indiceCount, 
        indicesVec.data(),
         GL_DYNAMIC_STORAGE_BIT
    );

    //  glCreateBuffers(1, &particleSystemDataBuffer);
    //  glNamedBufferStorage(
    //      particleSystemDataBuffer, 
    //      sizeof(ClothSystemDataBlock), 
    //      &particleSystemDataBlock,
    //       GL_DYNAMIC_STORAGE_BIT
    //  );


    glCreateBuffers(1, &edgesBuffer);
    glNamedBufferStorage(
       edgesBuffer, 
       sizeof(glm::ivec2) * edgeCount, 
       edgeVec.data(),
        GL_DYNAMIC_STORAGE_BIT
    );

    glCreateBuffers(1, &hingesBuffer);
    glNamedBufferStorage(
       hingesBuffer, 
       sizeof(glm::ivec4) * hingeCount, 
       hingeVec.data(),
        GL_DYNAMIC_STORAGE_BIT
    );

    glCreateBuffers(1, &edgeDebugDataBuffer);
    glNamedBufferStorage(
       edgeDebugDataBuffer, 
       sizeof(glm::dvec4) * 4, 
       edgeDebugDataVec,
        GL_DYNAMIC_STORAGE_BIT
    );

    glCreateBuffers(1, &edgeDebugEBO);
    glNamedBufferStorage(
       edgeDebugEBO, 
       sizeof(uint32_t) * 6, 
       edgeDebugEBOData,
        GL_DYNAMIC_STORAGE_BIT
    );
    
    glCreateBuffers(1, &pickedDebugDataBuffer);
    glNamedBufferStorage(
       pickedDebugDataBuffer, 
       sizeof(PickedClothData), 
       &pickedClothData,
        GL_DYNAMIC_STORAGE_BIT
    );

    glCreateBuffers(1, &fixedNodesDataBuffer);
    glNamedBufferStorage(
       fixedNodesDataBuffer, 
       sizeof(uint32_t) * fixedNodes.size(), 
       fixedNodes.data(),
        GL_DYNAMIC_STORAGE_BIT
    );

    glCreateBuffers(1, &clothTexcoordsBuffer);
    glNamedBufferStorage(
        clothTexcoordsBuffer,
        sizeof(glm::vec2) * clothTexcoords.size(),
        clothTexcoords.data(),
        GL_DYNAMIC_STORAGE_BIT
    );

    int texwidth, texheight, nrChannels;

    const std::string path = "./resources/assets/baka.jpg"; // Thanks C++.
    unsigned char *data = stbi_load(path.c_str(), &texwidth, &texheight, &nrChannels, 4);
    glCreateTextures(GL_TEXTURE_2D, 1, &clothTexture);

    glTextureParameteri(clothTexture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
    glTextureParameteri(clothTexture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
    glTextureParameteri(clothTexture, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTextureParameteri(clothTexture, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    glTextureStorage2D(clothTexture, 1, GL_RGBA8, texwidth, texheight);
    glTextureSubImage2D(clothTexture, 0, 0, 0, texwidth, texheight, GL_RGBA, GL_UNSIGNED_BYTE, data);
    stbi_image_free(data);

    //https://stackoverflow.com/questions/12399422/how-to-set-linker-flags-for-openmp-in-cmakes-try-compile-function

    glCreateVertexArrays(1, &fixedNodesVAO);
    glEnableVertexArrayAttrib(fixedNodesVAO, POSITION_ATTRIB_LOC);
    
    // nodeIdx attribute
    glVertexArrayAttribIFormat(fixedNodesVAO, POSITION_ATTRIB_LOC, 1, GL_UNSIGNED_INT, 0);
    glVertexArrayAttribBinding(fixedNodesVAO, POSITION_ATTRIB_LOC, 0);
    // Eigen::setNbThreads(12);


    #if PARDISO_SOLVE == 1
        // doPardiso();
        InitializePardiso();
    #endif
}

static Eigen::MatrixXd uvT(
    const Eigen::Ref< const Eigen::Vector3d> & u, 
    const Eigen::Ref< const Eigen::Vector3d> & v){
    return u * v.transpose();
}

static void hessTheta(
        const Eigen::Ref< const Eigen::Vector3d > & node0, 
        const Eigen::Ref< const Eigen::Vector3d > & node1, 
        const Eigen::Ref< const Eigen::Vector3d > & node2, 
        const Eigen::Ref< const Eigen::Vector3d > & node3, 
        double* hess){
    
    Eigen::Vector3d m_e0 = node1 - node0;
    Eigen::Vector3d m_e1 = node2 - node0;
    Eigen::Vector3d m_e2 = node3 - node0;
    Eigen::Vector3d m_e3 = node2 - node1;
    Eigen::Vector3d m_e4 = node3 - node1;

    double m_cosA1 = m_e0.dot(m_e1) / (m_e0.norm() * m_e1.norm());
    double m_cosA2 = m_e0.dot(m_e2) / (m_e0.norm() * m_e2.norm());
    double m_cosA3 = -m_e0.dot(m_e3) / (m_e0.norm() * m_e3.norm());
    double m_cosA4 = -m_e0.dot(m_e4) / (m_e0.norm() * m_e4.norm());

    double m_sinA1 = m_e0.cross(m_e1).norm() / (m_e0.norm() * m_e1.norm());
    double m_sinA2 = m_e0.cross(m_e2).norm() / (m_e0.norm() * m_e2.norm());
    double m_sinA3 = -m_e0.cross(m_e3).norm() / (m_e0.norm() * m_e3.norm());
    double m_sinA4 = -m_e0.cross(m_e4).norm() / (m_e0.norm() * m_e4.norm());

    Eigen::Vector3d m_nn1 = m_e0.cross(m_e3);
    m_nn1 = m_nn1 / m_nn1.norm();
    Eigen::Vector3d m_nn2 = -m_e0.cross(m_e4);
    m_nn2 = m_nn2 / m_nn2.norm();

    double m_h1 = m_e0.norm() * m_sinA1;
    double m_h2 = m_e0.norm() * m_sinA2;
    double m_h3 = -m_e0.norm() * m_sinA3;
    double m_h4 = -m_e0.norm() * m_sinA4;
    double m_h01 = m_e1.norm() * m_sinA1;
    double m_h02 = m_e2.norm() * m_sinA2;

    Eigen::Vector3d m_m1 = m_nn1.cross(m_e1) / m_e1.norm();
    Eigen::Vector3d m_m2 = -m_nn2.cross(m_e2) / m_e2.norm();
    Eigen::Vector3d m_m3 = -m_nn1.cross(m_e3) / m_e3.norm();
    Eigen::Vector3d m_m4 = m_nn2.cross(m_e4) / m_e4.norm();
    Eigen::Vector3d m_m01 = -m_nn1.cross(m_e0) / m_e0.norm();
    Eigen::Vector3d m_m02 = m_nn2.cross(m_e0) / m_e0.norm();
    //
    Eigen::Map<Eigen::MatrixXd> hess_theta(hess, 12, 12);
     hess_theta.setZero();

    Eigen::MatrixXd M331 = m_cosA3 / (m_h3 * m_h3) * uvT(m_m3, m_nn1);
    Eigen::MatrixXd M311 = m_cosA3 / (m_h3 * m_h1) * uvT(m_m1, m_nn1);
    Eigen::MatrixXd M131 = m_cosA1 / (m_h1 * m_h3) * uvT(m_m3, m_nn1);
    Eigen::MatrixXd M3011 = m_cosA3 / (m_h3 * m_h01) * uvT(m_m01, m_nn1);
    Eigen::MatrixXd M111 = m_cosA1 / (m_h1 * m_h1) * uvT(m_m1, m_nn1);
    Eigen::MatrixXd M1011 = m_cosA1 / (m_h1 * m_h01) * uvT(m_m01, m_nn1);
    
    Eigen::MatrixXd M442 = m_cosA4 / (m_h4 * m_h4) * uvT(m_m4, m_nn2);
    Eigen::MatrixXd M422 = m_cosA4 / (m_h4 * m_h2) * uvT(m_m2, m_nn2);
    Eigen::MatrixXd M242 = m_cosA2 / (m_h2 * m_h4) * uvT(m_m4, m_nn2);
    Eigen::MatrixXd M4022 = m_cosA4 / (m_h4 * m_h02) * uvT(m_m02, m_nn2);
    Eigen::MatrixXd M222 = m_cosA2 / (m_h2 * m_h2) * uvT(m_m2, m_nn2);
    Eigen::MatrixXd M2022 = m_cosA2 / (m_h2 * m_h02) * uvT(m_m02, m_nn2);

    Eigen::MatrixXd B1 = 1 / m_e0.norm() / m_e0.norm() * uvT(m_nn1, m_m01);
    Eigen::MatrixXd B2 = 1 / m_e0.norm() / m_e0.norm() * uvT(m_nn2, m_m02);

    Eigen::MatrixXd N13 = 1 / (m_h01 * m_h3) * uvT(m_nn1, m_m3);
    Eigen::MatrixXd N24 = 1 / (m_h02 * m_h4) * uvT(m_nn2, m_m4);
    Eigen::MatrixXd N11 = 1 / (m_h01 * m_h1) * uvT(m_nn1, m_m1);
    Eigen::MatrixXd N22 = 1 / (m_h02 * m_h2) * uvT(m_nn2, m_m2);
    Eigen::MatrixXd N101 = 1 / (m_h01 * m_h01) * uvT(m_nn1, m_m01);
    Eigen::MatrixXd N202 = 1 / (m_h02 * m_h02) * uvT(m_nn2, m_m02);

    hess_theta.block<3, 3>(0, 0) = mmt(M331) - B1 + mmt(M442) - B2;
    hess_theta.block<3, 3>(0, 3) = M311 + M131.transpose() + B1 + M422 + M242.transpose() + B2;
    hess_theta.block<3, 3>(0, 6) = M3011 - N13;
    hess_theta.block<3, 3>(0, 9) = M4022 - N24;
    hess_theta.block<3, 3>(3, 3) = mmt(M111) - B1 + mmt(M222) - B2;
    hess_theta.block<3, 3>(3, 6) = M1011 - N11;
    hess_theta.block<3, 3>(3, 9) = M2022 - N22;
    hess_theta.block<3, 3>(6, 6) = -mmt(N101);
    hess_theta.block<3, 3>(9, 9) = -mmt(N202);

    hess_theta.block<3, 3>(3, 0) = hess_theta.block<3, 3>(0, 3).transpose();
    hess_theta.block<3, 3>(6, 0) = hess_theta.block<3, 3>(0, 6).transpose();
    hess_theta.block<3, 3>(9, 0) = hess_theta.block<3, 3>(0, 9).transpose();
    hess_theta.block<3, 3>(6, 3) = hess_theta.block<3, 3>(3, 6).transpose();
    hess_theta.block<3, 3>(9, 3) = hess_theta.block<3, 3>(3, 9).transpose();
}

static void gradTheta(
        const Eigen::Ref< const Eigen::Vector3d >& node0, 
        const Eigen::Ref< const Eigen::Vector3d >& node1, 
        const Eigen::Ref< const Eigen::Vector3d >& node2, 
        const Eigen::Ref< const Eigen::Vector3d >& node3, 
        double* grad){
    
    Eigen::Vector3d m_e0 = node1 - node0;
    Eigen::Vector3d m_e1 = node2 - node0;
    Eigen::Vector3d m_e2 = node3 - node0;
    Eigen::Vector3d m_e3 = node2 - node1;
    Eigen::Vector3d m_e4 = node3 - node1;

    double m_cosA1 = m_e0.dot(m_e1) / (m_e0.norm() * m_e1.norm());
    double m_cosA2 = m_e0.dot(m_e2) / (m_e0.norm() * m_e2.norm());
    double m_cosA3 = -m_e0.dot(m_e3) / (m_e0.norm() * m_e3.norm());
    double m_cosA4 = -m_e0.dot(m_e4) / (m_e0.norm() * m_e4.norm());

    double m_sinA1 = m_e0.cross(m_e1).norm() / (m_e0.norm() * m_e1.norm());
    double m_sinA2 = m_e0.cross(m_e2).norm() / (m_e0.norm() * m_e2.norm());
    double m_sinA3 = -m_e0.cross(m_e3).norm() / (m_e0.norm() * m_e3.norm());
    double m_sinA4 = -m_e0.cross(m_e4).norm() / (m_e0.norm() * m_e4.norm());

    Eigen::Vector3d m_nn1 = m_e0.cross(m_e3);
    m_nn1 /= m_nn1.norm();
    Eigen::Vector3d m_nn2 = -m_e0.cross(m_e4);
    m_nn2 /= m_nn2.norm();

    double m_h1 = m_e0.norm() * m_sinA1;
    double m_h2 = m_e0.norm() * m_sinA2;
    double m_h3 = -m_e0.norm() * m_sinA3;
    double m_h4 = -m_e0.norm() * m_sinA4;
    double m_h01 = m_e1.norm() * m_sinA1;
    double m_h02 = m_e2.norm() * m_sinA2;

    Eigen::Vector3d gradTheta0 = m_cosA3 * m_nn1 / m_h3 + m_cosA4 * m_nn2 / m_h4;
    Eigen::Vector3d gradTheta1 = m_cosA1 * m_nn1 / m_h1 + m_cosA2 * m_nn2 / m_h2;
    Eigen::Vector3d gradTheta2 = -m_nn1 / m_h01;
    Eigen::Vector3d gradTheta3 = -m_nn2 / m_h02;

    Eigen::Map<Eigen::Matrix<double, 12, 1>> dFMap(grad);
    dFMap.segment<3>(0) = gradTheta0;
    dFMap.segment<3>(3) = gradTheta1;
    dFMap.segment<3>(6) = gradTheta2;
    dFMap.segment<3>(9) = gradTheta3;
}

static void subFromSparseEntries(const Eigen::Ref< const Eigen::MatrixXd > & mat, 
                                uint32_t* indices, 
                                SparseEntries* sparseEntries,
                                const std::span<uint32_t> oldToNewIndices){
    for (uint32_t i = 0; i < mat.rows(); ++i){
        uint32_t oldRowIdx = indices[i];
        uint32_t newRowIdx = oldToNewIndices[oldRowIdx];
        if (newRowIdx == UINT32_MAX) {
            continue;
        }

        for (uint32_t j = 0; j < mat.cols(); ++j){
            uint32_t oldColIdx = indices[j];
            uint32_t newColIdx = oldToNewIndices[oldColIdx];
            if (newColIdx == UINT32_MAX) {
                continue;
            }
            
            #if PARDISO_SOLVE == 1
            if (newColIdx < newRowIdx) {
                continue;
            }
            #endif

            if (mat(i, j) == 0.0) {
                continue;
            }
            
            sparseEntries->push_back(Eigen::Triplet<double>(newRowIdx, newColIdx, -mat(i, j)));
        }
    }
}

static void calculateGradEb_HessEb_Shell(
    const Eigen::Ref< const Eigen::Vector3d > & node0, 
    const Eigen::Ref< const Eigen::Vector3d > & node1, 
    const Eigen::Ref< const Eigen::Vector3d > & node2, 
    const Eigen::Ref< const Eigen::Vector3d > & node3, 
    double thetaBar, double kb, double* dF, double* dJ){

    double theta = getTheta(node0, node1, node2, node3);
    double grad[12];
    gradTheta(node0, node1, node2, node3, grad);
    
    Eigen::Map<
        Eigen::Matrix<
            double,
            12, 1>
    > gradMap(grad);
    Eigen::Map<
        Eigen::Matrix<
            double,
            12, 1>
    > dFMap(dF);

    dFMap.segment<12>(0) = 0.5 * kb * (2.0 * (theta - thetaBar) * gradMap);
    
    double hess[144];
    hessTheta(node0, node1, node2, node3, hess);
    
    
    Eigen::Map<Eigen::Matrix<double, 12, 12>> dJMap(dJ);
    Eigen::Map<Eigen::Matrix<double, 12, 12>> hessMap(hess);

    Eigen::MatrixXd gradgradT = gradMap*gradMap.transpose() ;
    dJMap.block<12, 12>(0, 0) = 0.5 * kb * (2.0 * gradgradT + 2.0 * (theta - thetaBar) * hessMap);
}

static void calculateFb_Jb(double* bendingForcesVecPtr, SparseEntries* sparseEntries, // double* bendingHessiansVecPtr, 
                        const std::vector<glm::ivec4>& hingeVec, 
                        const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 1>>& q, 
                        double kb, 
                        const std::vector<double>& thetaBar, 
                        uint32_t dofCount, uint32_t hingeCount,
                        const std::span<uint32_t> oldToNewIndiceMapping) {
    
    Eigen::Map<Eigen::Matrix<
            double,
            Eigen::Dynamic, 1
            >> bendingForces(bendingForcesVecPtr, dofCount, 1);

    // Bending Forces
    for (uint32_t hingeIdx = 0; hingeIdx < hingeCount; ++hingeIdx) {
        glm::ivec4 kHinge = hingeVec[hingeIdx];
        uint32_t node0 = kHinge[0];
        uint32_t node1 = kHinge[1];
        uint32_t node2 = kHinge[2];
        uint32_t node3 = kHinge[3];

        Eigen::Vector3d x0 = q.segment<3>(3 * node0);
        Eigen::Vector3d x1 = q.segment<3>(3 * node1);
        Eigen::Vector3d x2 = q.segment<3>(3 * node2);
        Eigen::Vector3d x3 = q.segment<3>(3 * node3);

        uint32_t ind[] = { 3 * node0, 3 * node0 + 1, 3 * node0 + 2,
                      3 * node1, 3 * node1 + 1, 3 * node1 + 2,
                      3 * node2, 3 * node2 + 1, 3 * node2 + 2,
                      3 * node3, 3 * node3 + 1, 3 * node3 + 2 };

        double dF[12];
        double dJ[144];


        Eigen::Map<
        Eigen::Matrix<
            double,
            12, 1
        >> dFMap(dF);

        Eigen::Map<
            Eigen::Matrix<
            double,
            12, 12>
        > dJMap(dJ);
        
        calculateGradEb_HessEb_Shell(x0, x1, x2, x3, thetaBar[hingeIdx], kb, dF, dJ);
        bendingForces(ind) = bendingForces(ind).eval() - dFMap;
        subFromSparseEntries(dJMap, ind, sparseEntries, oldToNewIndiceMapping);
    }
}


static void calculateGradEs_HessEs(
    const Eigen::Ref< const Eigen::Vector3d>& node0, 
    const Eigen::Ref< const Eigen::Vector3d >& node1,
    double lk, double EA, double* dF, double* dJ) {

   // Gradient of Es
    Eigen::Vector3d edge = node1 - node0;
    double edgeLen = edge.norm();
    Eigen::Vector3d tangent = edge / edgeLen;
    double epsX = edgeLen / lk - 1.0;
    Eigen::Vector3d dF_unit = EA * tangent * epsX;
    
    Eigen::Map<
        Eigen::Matrix<
            double, 
            6, 1>
    > dFMap(dF);
    dFMap << -dF_unit, dF_unit;

    Eigen::Matrix3d Id3 = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d M = EA * ((1 / lk - 1 / edgeLen) * Id3 + 1 / edgeLen * (edge * edge.transpose()) / (edgeLen * edgeLen));
    
    Eigen::Map<Eigen::Matrix<double, 6, 6>> dJMap(dJ);
    dJMap.setZero();
    dJMap.block<3, 3>(0, 0) = M;
    dJMap.block<3, 3>(3, 3) = M;
    dJMap.block<3, 3>(0, 3) = -M;
    dJMap.block<3, 3>(3, 0) = -M;    
}

static void calculateFs_Js(double* stretchingForcesVecPtr, SparseEntries* sparseEntries,// double* stretchingHessiansVecPtr, 
            const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 1>>& q, 
            const std::vector<glm::ivec2>& edgeVec, 
            const std::vector<double>& undeformedEdgeLengthVec, 
            const std::vector<double>& elasticStretchingVec, 
            uint32_t dofCount, uint32_t edgeCount,
            const std::span<uint32_t> oldToNewIndiceMapping) {

    Eigen::Map<
        Eigen::Matrix<
            double,
            Eigen::Dynamic, 1>
    > stretchForces(stretchingForcesVecPtr, dofCount, 1);

    for (uint32_t edgeIdx = 0; edgeIdx < edgeCount; edgeIdx++) {
        glm::ivec2 kEdge = edgeVec[edgeIdx];
        uint32_t node0 = kEdge[0];
        uint32_t node1 = kEdge[1];
        double lk = undeformedEdgeLengthVec[edgeIdx];
        double ks  = elasticStretchingVec[edgeIdx];

        Eigen::Vector3d x0 = q.segment<3>(3 * node0);
        Eigen::Vector3d x1 = q.segment<3>(3 * node1);

        uint32_t ind[] = { 3 * node0, 3 * node0 + 1, 3 * node0 + 2,
                      3 * node1, 3 * node1 + 1, 3 * node1 + 2 };
        
        double dF[6];
        double dJ[36];
        Eigen::Map<
            Eigen::Matrix<
                double, 
                6, 1>
        > dFMap(dF);
        Eigen::Map<
            Eigen::Matrix<
                double, 
                6, 6>
        > dJMap(dJ);

        calculateGradEs_HessEs(x0, x1, lk, ks, dF, dJ);
        stretchForces(ind) = stretchForces(ind).eval() - dFMap;
        subFromSparseEntries(dJMap, ind, sparseEntries, oldToNewIndiceMapping);
    }
}

static void SolveSystemWithConjugateGradient(const Eigen::Ref<const Eigen::SparseMatrix<double>>& J_free,
                        const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 1>>& f_free,
                        Eigen::Ref<Eigen::Matrix<double, Eigen::Dynamic, 1>> dq_free) {
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
    solver.compute(J_free);
    dq_free = solver.solve(f_free);
}

#if PARDISO_SOLVE == 1

static void shiftIndices(int n, int nonzeros, int* ia, int* ja, int value)
{
  int i;
    for (i = 0; i < n+1; i++) 
    {
        ia[i] += value;
    }
    for (i = 0; i < nonzeros; i++) 
    {
        ja[i] += value;
    }
}

void ClothSystem::InitializePardiso(){
    pardisoData.mtype = -2;        /* Real symmetric */
    /* -------------------------------------------------------------------- */
    /* ..  Setup Pardiso control parameters and initialize the solvers      */
    /*     internal adress pointers. This is only necessary for the FIRST   */
    /*     call of the PARDISO solver.                                      */
    /* ---------------------------------------------------------------------*/

    pardisoData.error = 0;
    pardisoData.solver = 0; /* use sparse direct solver */
    pardisoinit_d(pardisoData.pt, &pardisoData.mtype, &pardisoData.solver, &pardisoData.iparm[1], pardisoData.dparm, &pardisoData.error);

    pardisoData.num_procs = 12;
    pardisoData.iparm[3] = pardisoData.num_procs;
    pardisoData.iparm[7] = 1;

    pardisoData.maxfct = 1;         /* Maximum number of numerical factorizations.  */
    pardisoData.mnum = 1;         /* Which factorization to use. */

    pardisoData.msglvl = 0;         /* Print statistical information  */
    pardisoData.error = 0;         /* Initialize error flag */
            
    pardisoData.n = freeDOFCount;         /* Number of equations */
    pardisoData.nrhs = 1;         /* Number of right hand sides. */
}


static void SolveSystemWithPardiso(Eigen::Ref< Eigen::SparseMatrix<double, Eigen::RowMajor>> J_free,
    Eigen::Ref<Eigen::Matrix<double, Eigen::Dynamic, 1>> f_free,
    Eigen::Ref<Eigen::Matrix<double, Eigen::Dynamic, 1>> dq_free,
    ClothSystem::PardisoData& pardisoData) {

    int    *ia = J_free.outerIndexPtr();
    int    *ja = J_free.innerIndexPtr();

    double  *a = J_free.valuePtr();

    double  *b = f_free.data();
    double  *x = dq_free.data();

    int      nnz = ia[pardisoData.n];

    /* -------------------------------------------------------------------- */
    /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
    /*     notation.                                                        */
    /* -------------------------------------------------------------------- */
    shiftIndices(pardisoData.n, nnz, ia, ja, 1);


    /* -------------------------------------------------------------------- */
    /*  .. pardiso_chk_matrix(...)                                          */
    /*     Checks the consistency of the given matrix.                      */
    /*     Use this functionality only for debugging purposes               */
    /* -------------------------------------------------------------------- */
    // pardiso_chkmatrix_d(&pardisoData.mtype, &pardisoData.n, a, ia, ja, &pardisoData.error);
    // if (pardisoData.error != 0)
    // {
    //     std::printf("\nERROR in consistency of matrix: %d", pardisoData.error);
    // }


    // /* -------------------------------------------------------------------- */    
    // /* ..  Reordering and Symbolic Factorization.  This step also allocates */
    // /*     all memory that is necessary for the factorization.              */
    // /* -------------------------------------------------------------------- */ 
    
    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */
    pardisoData.phase = 13;


    pardiso_d(pardisoData.pt, &pardisoData.maxfct, &pardisoData.mnum, &pardisoData.mtype, &pardisoData.phase,
        &pardisoData.n, a, ia, ja, &idum, &pardisoData.nrhs,
        &pardisoData.iparm[1], &pardisoData.msglvl, b, x, &pardisoData.error, pardisoData.dparm);

    if (pardisoData.error != 0)
    {
        std::printf("\nERROR during solve: %d", pardisoData.error);
        return;
    }

    /* -------------------------------------------------------------------- */
    /* ..  Convert matrix from 1-based Fortan notation to 0-based C         */
    /*     notation.                                                        */
    /* -------------------------------------------------------------------- */
    shiftIndices(pardisoData.n, nnz, ia, ja, -1);

    /* -------------------------------------------------------------------- */
    /* ..  Termination and release of memory.                               */
    /* -------------------------------------------------------------------- */
    pardisoData.phase = -1;                 /* Release internal memory. */

    pardiso_d(pardisoData.pt, &pardisoData.maxfct, &pardisoData.mnum, &pardisoData.mtype, &pardisoData.phase,
        &pardisoData.n, a, ia, ja, &idum, &pardisoData.nrhs,
        &pardisoData.iparm[1], &pardisoData.msglvl, &ddum, &ddum, &pardisoData.error, pardisoData.dparm);
}

#endif

void ClothSystem::CalculateForces() {

    auto currTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = currTime - lastTime;

    
    if (duration.count() > dt ) {
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> bufferStrides(4, 1);
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> tmpVecStrides(3, 1);
        Eigen::Map<
            Eigen::Matrix<
                double, 
                Eigen::Dynamic, 1>
        > massVec(massVector.data(), dofCount, 1);

        uint32_t iter = 0;
        double error = 10 * tol;

        Eigen::Map<
            Eigen::Matrix<
                double, 
                Eigen::Dynamic, 1>
        > externalForces(externalForcesVec.data(), dofCount, 1);

        Eigen::Map<
            Eigen::Matrix<
                double, 
                Eigen::Dynamic, 1>
        > q0(lastdofPositions.data(), dofCount, 1);
        Eigen::Map<
            Eigen::Matrix<
                double, 
                Eigen::Dynamic, 1>
        > q(dofPositions.data(), dofCount, 1);
        Eigen::Map<
            Eigen::Matrix<
                double, 
                Eigen::Dynamic, 1>
        > u(dofVelocities.data(), dofCount, 1);

        q0 = q.eval();

        Eigen::Matrix<double, Eigen::Dynamic, 1> dq_free(freeDOFCount);

        while (error > tol) {
            if (iter > 5) { // prevent lagging
                break;
            }
            
            auto startForceTime = std::chrono::high_resolution_clock::now();
            std::vector<double> bendingForcesVec(dofCount, 0.0);

            std::vector<double> stretchingForcesVec(dofCount, 0.0);

            Eigen::Map<
                Eigen::Matrix<
                    double, 
                    Eigen::Dynamic, 1>
            > bendingForces(bendingForcesVec.data(), dofCount, 1);

            Eigen::Map<
                Eigen::Matrix<
                    double, 
                    Eigen::Dynamic, 1>
            > stretchForces(stretchingForcesVec.data(), dofCount, 1);

            Eigen::SparseMatrix<double, 
            #if PARDISO_SOLVE == 1
             Eigen::RowMajor
            #else
             Eigen::ColMajor
            #endif
            > jacobianForces(freeDOFCount, freeDOFCount);
            SparseEntries hessianEntries;
            
            // Bending Forces -------------------------------------------------------------
            calculateFb_Jb(bendingForcesVec.data(), &hessianEntries, hingeVec, q, bendingStiffness, thetaBarVec, dofCount, hingeCount, oldToNewIndiceMapping);
            
            // Stretching Forces -------------------------------------------------------------
            calculateFs_Js(stretchingForcesVec.data(), &hessianEntries, q, edgeVec, undeformedEdgeLengthVec, elasticStretchingVec, dofCount, edgeCount, oldToNewIndiceMapping);
           
            auto endForceTime = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> forceTime = endForceTime - startForceTime;
            totalForceTime += forceTime.count();
            
            auto startSolverTime = std::chrono::high_resolution_clock::now();

            jacobianForces.setFromTriplets(hessianEntries.begin(), hessianEntries.end());
            Eigen::Matrix<double, Eigen::Dynamic, 1> totalForces = bendingForces + stretchForces + externalForces;
            Eigen::Matrix<double, Eigen::Dynamic, 1> f = (massVec / dt).cwiseProduct((q - q0) / dt - u) - totalForces;
            Eigen::Matrix<double, Eigen::Dynamic, 1> f_free = f(freeIndices);

            Eigen::SparseMatrix<double,
            #if PARDISO_SOLVE == 1
             Eigen::RowMajor
            #else
             Eigen::ColMajor
            #endif
             > J_free = massMatrix  - jacobianForces;
            J_free.makeCompressed();
            
            #if PARDISO_SOLVE == 1
            SolveSystemWithPardiso(J_free, f_free, dq_free, pardisoData);
            #else
            SolveSystemWithConjugateGradient(J_free, f_free, dq_free);
            #endif
            auto endSolverTime = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> solverTime = endSolverTime - startSolverTime;      
            totalSolveTime += solverTime.count();  
            
            q(freeIndices) = q(freeIndices).eval() - dq_free;

            error = f_free.cwiseAbs().sum();
            //error = 0.0;
            iter += 1;
        }

        //fmt::print("Time: {:.6f}\n[", totalTime);
        //uint32_t perline = 0;
        //for (uint32_t freeIdx : freeIndices) {
        //   fmt::print("{: .8e} ", q(freeIdx));
        //   perline += 1;
        //   if (perline == 4) {
        //       perline = 0;
        //       fmt::print("\n ");
        //   }
        //}
        //fmt::print("]\n");

        u = (q - q0) / dt;

        for (uint32_t particleIdx = 0; particleIdx < particleCount; particleIdx++) {
            clothPositionVec[particleIdx * 4    ] = dofPositions[particleIdx * 3    ];
            clothPositionVec[particleIdx * 4 + 1] = dofPositions[particleIdx * 3 + 1];
            clothPositionVec[particleIdx * 4 + 2] = dofPositions[particleIdx * 3 + 2];
        }

        glNamedBufferSubData(positionBuffer, 0, sizeof(float) * 4 * particleCount, clothPositionVec.data());
        lastTime = currTime;
        totalTime += dt;
        
    }
}

void ClothSystem::BindBuffers() const{
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_POSITIONS_SSBO_BINDING, positionBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_SYSTEM_DATA_SSBO_BINDING, particleSystemDataBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_EDGE_DATA_SSBO_BINDING, edgesBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_HINGE_DATA_SSBO_BINDING, hingesBuffer);
    glBindBufferBase(GL_UNIFORM_BUFFER, PICKED_DATA_UBO_BINDING, pickedDebugDataBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_FIXED_NODES_SSBO_BINDING, fixedNodesDataBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_TEXCOORDS_SSBO_BINDING, clothTexcoordsBuffer);
    glBindTexture(GL_TEXTURE_2D, clothTexture);

    particleComputeShader->use();
}

void ClothSystem::SetupRender() {

    particleShader->use();
}


void ClothSystem::RenderDebug(GLuint VAO) {
    // glBindVertexArray(VAO);
    // glNamedBufferSubData(pickedDebugDataBuffer, 0, sizeof(PickedClothData), &pickedClothData);
    glNamedBufferSubData(fixedNodesDataBuffer, 0, sizeof(uint32_t)* fixedNodesCount, fixedNodes.data());

    //edgeDebugShader->use();
    //glVertexArrayVertexBuffer(VAO, 0, edgeDebugDataBuffer, 0, sizeof(glm::dvec4));
    //glVertexArrayElementBuffer(VAO, edgeDebugEBO);
    //glPointSize(10.0);
    //glLineWidth(10.0);
    //glDrawElements(GL_LINES, 2, GL_UNSIGNED_INT, 0);
    //glDrawElements(GL_POINTS, 2, GL_UNSIGNED_INT, 0);

    //hingeDebugShader->use();
    //glPointSize(20.0);
    //glLineWidth(20.0);
    //glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
    //glDrawElements(GL_POINTS, 6, GL_UNSIGNED_INT, 0);
    
    
    glBindVertexArray(fixedNodesVAO);
    fixedNodesDebugShader->use();
    glPointSize(20.0);
    glVertexArrayVertexBuffer(fixedNodesVAO, 0, fixedNodesDataBuffer, 0, sizeof(uint32_t));
    glDrawArrays(GL_POINTS, 0, fixedNodesCount);


    auto currTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = currTime - lastTime;
    
    //if (duration.count() > 120) {
    //    pickedClothData.pickedEdgeIdx = (pickedClothData.pickedEdgeIdx + 1) % edgeCount;
    //    // fmt::print("Edge {}: {} {}\n",
    //    //     pickedClothData.pickedEdgeIdx, 
    //    //     edgeVec[pickedClothData.pickedEdgeIdx].x, 
    //    //     edgeVec[pickedClothData.pickedEdgeIdx].y);

    //    // fmt::print("Hinge {}: {} {} {} {}\n",
    //    //     pickedClothData.pickedHingeIdx, 
    //    //     hingeVec[pickedClothData.pickedHingeIdx].x, 
    //    //     hingeVec[pickedClothData.pickedHingeIdx].y,
    //    //     hingeVec[pickedClothData.pickedHingeIdx].z,
    //    //     hingeVec[pickedClothData.pickedHingeIdx].w);
    //    pickedClothData.pickedHingeIdx = (pickedClothData.pickedHingeIdx + 1) % hingeCount;        
    //    lastTime = currTime;
    //}
}