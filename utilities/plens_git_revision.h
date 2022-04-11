#ifndef PLENS_GIT_REVISION_H
#define PLENS_GIT_REVISION_H

// This file will be used by `plens/CMakeLists.txt` (the main cmake script) to generate
// `plens_git_revision.h` in this directory. The latter file can be used to print git-related info
// from the code.

/**
 * Name of the local git branch of the source directory.
 */
#define PLENS_GIT_BRANCH "entropy_var_grad"

/**
 * Full sha1 revision of the current git HEAD.
 */
#define PLENS_GIT_REVISION "f9846c1118d7e60a48db6fb9d4d48b9e8273cc20"

/**
 * Short sha1 revision of the current git HEAD.
 */
#define PLENS_GIT_SHORTREV "f9846c1"

#endif
