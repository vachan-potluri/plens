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
#define PLENS_GIT_REVISION "e804bb708707f0ebd2281d50668cc455f22176d4"

/**
 * Short sha1 revision of the current git HEAD.
 */
#define PLENS_GIT_SHORTREV "e804bb7"

#endif
