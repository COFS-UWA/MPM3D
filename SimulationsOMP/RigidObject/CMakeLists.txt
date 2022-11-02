target_sources(
    SimulationsOMP PRIVATE
    # 2D contact
    ContactModel2D.h ContactModel2D.cpp
    NoContact2D.h NoContact2D.cpp
    SmoothContact2D.h SmoothContact2D.cpp
    FrictionalContact2D.h FrictionalContact2D.cpp
    StickyContact2D.h StickyContact2D.cpp
    RoughContact2D.h RoughContact2D.cpp
    Force2D.h
    RigidCircle.h RigidCircle.cpp
    RigidObjectMotion2D.h RigidObjectMotion2D.cpp
    RigidMeshT2D.h RigidMeshT2D.cpp
    RigidObjectByT2DMesh.h RigidObjectByT2DMesh.cpp
    # 3D contact
    ContactModel3D.h ContactModel3D.cpp
    NoContact3D.h NoContact3D.cpp
    SmoothContact3D.h SmoothContact3D.cpp
    SmoothContact3D_Quad.h SmoothContact3D_Quad.cpp
    RoughContact3D.h RoughContact3D.cpp
    FrictionalContact3D.h FrictionalContact3D.cpp
    StickyContact3D.h StickyContact3D.cpp
	NonSeperateContact3D.h NonSeperateContact3D.cpp
    Force3D.h
    RigidCube.h RigidCube.cpp
    RigidCylinder.h RigidCylinder.cpp
    RigidCone.h RigidCone.cpp
    #
    RigidObjectMotion3D.h RigidObjectMotion3D.cpp
    RigidMeshT3D.h RigidMeshT3D.cpp
    RigidObjectByT3DMesh.h RigidObjectByT3DMesh.cpp
    # utilities
    RigidObject_hdf5_utilities.h RigidObject_hdf5_utilities.cpp
    )