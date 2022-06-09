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
#define PLENS_GIT_REVISION "444c33fd0424b5c4352630b675f520bb726eb3b2"

/**
 * Short sha1 revision of the current git HEAD.
 */
#define PLENS_GIT_SHORTREV "444c33f"

#endif
