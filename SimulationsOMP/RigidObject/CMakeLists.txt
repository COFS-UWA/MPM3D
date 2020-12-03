target_sources(
    SimulationsOMP PRIVATE
    ContactModel3D.h ContactModel3D.cpp
    SmoothContact3D.h SmoothContact3D.cpp
    RoughContact3D.h RoughContact3D.cpp
    FrictionalContact3D.h FrictionalContact3D.cpp
    Force2D.h Force3D.h
    RigidCube.h RigidCube.cpp
    RigidCylinder.h RigidCylinder.cpp
    RigidCone.h RigidCone.cpp
    RigidObject_hdf5_utilities.h RigidObject_hdf5_utilities.cpp
    )