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
#define PLENS_GIT_REVISION "e5aea21c7a1f1e52acc2678b6c7340b78368a0cb"

/**
 * Short sha1 revision of the current git HEAD.
 */
#define PLENS_GIT_SHORTREV "e5aea21"

#endif
