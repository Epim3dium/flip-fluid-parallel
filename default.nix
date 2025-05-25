with import <nixpkgs> {};

mkShell {
    name = "vulkan";
    packages = [
        glm
        glfw
        freetype
        valgrind
        kdePackages.kcachegrind
    ];

    buildInputs = with pkgs; [
        glm
        glfw
        sfml
        freetype
        vulkan-loader
    ];

    SFML_PATH = "${sfml}/lib/cmake";
    shellHook = ''
        # Prepend or append directories to the PATH
        export PATH="$PATH:$LD_LIBRARY_PATH:$SFML_PATH"
    '';
}
