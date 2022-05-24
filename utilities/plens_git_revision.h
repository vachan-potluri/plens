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
#define PLENS_GIT_REVISION "2bff91cc8708c33e50a2e5d6dc9ef7d59947e62a"

/**
 * Short sha1 revision of the current git HEAD.
 */
#define PLENS_GIT_SHORTREV "2bff91c"

#endif
