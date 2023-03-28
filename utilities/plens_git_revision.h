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
#define PLENS_GIT_REVISION "52905eb904bd0f0a065b77a2ecf3dff78d777a6b"

/**
 * Short sha1 revision of the current git HEAD.
 */
#define PLENS_GIT_SHORTREV "52905eb"

#endif
