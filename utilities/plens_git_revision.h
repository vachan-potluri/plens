#ifndef PLENS_GIT_REVISION_H
#define PLENS_GIT_REVISION_H

// This file will be used by `plens/CMakeLists.txt` (the main cmake script) to generate
// `plens_git_revision.h` in this directory. The latter file can be used to print git-related info
// from the code.

/**
 * Name of the local git branch of the source directory.
 */
#define PLENS_GIT_BRANCH "dif_residual_correction"

/**
 * Full sha1 revision of the current git HEAD.
 */
#define PLENS_GIT_REVISION "c8c91a629cf9c7a54b5e61cde45a0e6648c653b5"

/**
 * Short sha1 revision of the current git HEAD.
 */
#define PLENS_GIT_SHORTREV "c8c91a6"

#endif
