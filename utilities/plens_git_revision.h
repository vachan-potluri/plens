#ifndef PLENS_GIT_REVISION_H
#define PLENS_GIT_REVISION_H

// This file will be used by `plens/CMakeLists.txt` (the main cmake script) to generate
// `plens_git_revision.h` in this directory. The latter file can be used to print git-related info
// from the code.

/**
 * Name of the local git branch of the source directory.
 */
#define PLENS_GIT_BRANCH "filtering"

/**
 * Full sha1 revision of the current git HEAD.
 */
#define PLENS_GIT_REVISION "cf1a6112da2d7486d07c1c46e1d4d949ee127d43"

/**
 * Short sha1 revision of the current git HEAD.
 */
#define PLENS_GIT_SHORTREV "cf1a611"

#endif
