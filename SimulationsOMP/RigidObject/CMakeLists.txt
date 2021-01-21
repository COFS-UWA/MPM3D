target_sources(
    SimulationsOMP PRIVATE
    # 2D contact
    Force2D.h
    RigidCircle.h RigidCircle.cpp
    # 3D contact
    ContactModel3D.h ContactModel3D.cpp
    NoContact3D.h NoContact3D.cpp
    SmoothContact3D.h SmoothContact3D.cpp
    RoughContact3D.h RoughContact3D.cpp
    FrictionalContact3D.h FrictionalContact3D.cpp
    Force3D.h
    RigidCube.h RigidCube.cpp
    RigidCylinder.h RigidCylinder.cpp
    RigidCone.h RigidCone.cpp
    RigidObjectMotion3D.h RigidObjectMotion3D.cpp
    RigidMeshT3D.h RigidMeshT3D.cpp
    RigidObjectByT3DMesh.h RigidObjectByT3DMesh.cpp
    # utilities
    RigidObject_hdf5_utilities.h RigidObject_hdf5_utilities.cpp
    )