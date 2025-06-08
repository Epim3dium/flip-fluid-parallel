with import <nixpkgs> {};

let
    gcc12 = pkgs.gcc12;
    sfmlWithGcc12 = pkgs.sfml.overrideAttrs (old: {
        nativeBuildInputs = old.nativeBuildInputs ++ [ pkgs.gcc12 ];
        CXX = "${pkgs.gcc12}/bin/g++";
        CC  = "${pkgs.gcc12}/bin/gcc";
    });
in
mkShell {
    name = "cuda-and-openmp";
    stdenv = gcc12Stdenv;
    packages = [
        glm
        glfw
        freetype
        valgrind
        kdePackages.kcachegrind
    ];

    buildInputs = with pkgs; [
        cudaPackages.cuda_cudart
        cudaPackages.cuda_nvcc
        cudaPackages.cuda_cccl
        cudaPackages.cudatoolkit
        cudaPackages.cuda_nvprof
        cudaPackages.cuda_nvvp
        jdk8
        linuxPackages.nvidia_x11
        gcc12
        libGLU libGL
        glm
        glfw
        sfmlWithGcc12
        freetype
        vulkan-loader
    ];

    SFML_PATH = "${sfmlWithGcc12}/lib/cmake";
    shellHook = ''
        export JAVA_HOME=${pkgs.jdk8}
        export CC=${pkgs.gcc12}/bin/gcc
        export CXX=${pkgs.gcc12}/bin/g++

        export CUDAHOSTCXX=${pkgs.gcc12}/bin/g++
        export CUDA_HOST_COMPILER=${pkgs.gcc12}/bin/gcc

        export CUDA_HOME=${pkgs.cudaPackages.cuda_cudart}
        export CUDA_PATH=${pkgs.cudaPackages.cuda_cudart}

        export LD_LIBRARY_PATH=${pkgs.cudaPackages.cuda_cudart}/lib64:${pkgs.cudaPackages.cuda_cudart}/lib:$LD_LIBRARY_PATH
        export LD_LIBRARY_PATH=${stdenv.cc.cc.lib}/lib:$LD_LIBRARY_PATH
        export LD_LIBRARY_PATH=${pkgs.linuxPackages.nvidia_x11}/lib:$LD_LIBRARY_PATH
        export LD_LIBRARY_PATH="$JAVA_HOME/lib:$JAVA_HOME/lib/server:$LD_LIBRARY_PATH"

        export LIBRARY_PATH=${pkgs.cudaPackages.cuda_cudart}/lib64:${pkgs.cudaPackages.cuda_cudart}/lib:$LIBRARY_PATH
            export PATH=$SFML_PATH/bin:$PATH
    '';
}
