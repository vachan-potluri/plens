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
#define PLENS_GIT_REVISION "9fac1c6964cb7ad8ad249724fa0a0eeb3c2fc80b"

/**
 * Short sha1 revision of the current git HEAD.
 */
#define PLENS_GIT_SHORTREV "9fac1c6"

#endif
