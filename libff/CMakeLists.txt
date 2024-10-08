set(FF_EXTRASRCS)
set(FF_EXTRALIBS)
if(${CURVE} STREQUAL "BN128")
  include_directories(depends)
  set(
    FF_EXTRASRCS

    ${FF_EXTRASRCS}
    algebra/curves/bn128/bn128_fields.cpp
    algebra/curves/bn128/bn128_g1.cpp
    algebra/curves/bn128/bn128_g2.cpp
    algebra/curves/bn128/bn128_gt.cpp
    algebra/curves/bn128/bn128_init.cpp
    algebra/curves/bn128/bn128_pairing.cpp
    algebra/curves/bn128/bn128_pp.cpp
  )
  set(
    FF_EXTRALIBS

    ${FF_EXTRALIBS}
    zm
  )
endif()


set(FF_BINARY_FIELDS)
if("${BINARY_FIELDS}")
  add_definitions(
          -DBINARY_FIELDS
  )

  set(FF_BINARY_FIELDS

     ${FF_BINARY_FIELDS}
      algebra/fields/binary/gf32.cpp
      algebra/fields/binary/gf64.cpp
      algebra/fields/binary/gf128.cpp
      algebra/fields/binary/gf192.cpp
      algebra/fields/binary/gf256.cpp
  )
endif()

add_library(
  ff
  STATIC

  algebra/curves/bls12_381/bls12_381_fields.cpp
  algebra/curves/bls12_381/bls12_381_g1.cpp
  algebra/curves/bls12_381/bls12_381_g2.cpp
  algebra/curves/bls12_381/bls12_381_init.cpp
  algebra/curves/bls12_381/bls12_381_pairing.cpp
  algebra/curves/bls12_381/bls12_381_pp.cpp

  algebra/curves/alt_bn128/alt_bn128_fields.cpp
  algebra/curves/alt_bn128/alt_bn128_g1.cpp
  algebra/curves/alt_bn128/alt_bn128_g2.cpp
  algebra/curves/alt_bn128/alt_bn128_init.cpp
  algebra/curves/alt_bn128/alt_bn128_pairing.cpp
  algebra/curves/alt_bn128/alt_bn128_pp.cpp

  algebra/curves/bn124/bn124_fields.cpp
  algebra/curves/bn124/bn124_g1.cpp
  algebra/curves/bn124/bn124_g2.cpp
  algebra/curves/bn124/bn124_init.cpp
  algebra/curves/bn124/bn124_pairing.cpp
  algebra/curves/bn124/bn124_pp.cpp

  algebra/curves/bn183/bn183_fields.cpp
  algebra/curves/bn183/bn183_g1.cpp
  algebra/curves/bn183/bn183_g2.cpp
  algebra/curves/bn183/bn183_init.cpp
  algebra/curves/bn183/bn183_pairing.cpp
  algebra/curves/bn183/bn183_pp.cpp

  algebra/curves/bn254/bn254_fields.cpp
  algebra/curves/bn254/bn254_g1.cpp
  algebra/curves/bn254/bn254_g2.cpp
  algebra/curves/bn254/bn254_init.cpp
  algebra/curves/bn254/bn254_pairing.cpp
  algebra/curves/bn254/bn254_pp.cpp

  algebra/curves/baby_jubjub/baby_jubjub_fields.cpp
  algebra/curves/baby_jubjub/baby_jubjub_g1.cpp
  algebra/curves/baby_jubjub/baby_jubjub_init.cpp
  algebra/curves/baby_jubjub/baby_jubjub_pp.cpp

  algebra/curves/jubjub_bn124/jubjub_bn124_fields.cpp
  algebra/curves/jubjub_bn124/jubjub_bn124_g1.cpp
  algebra/curves/jubjub_bn124/jubjub_bn124_init.cpp
  algebra/curves/jubjub_bn124/jubjub_bn124_pp.cpp

  algebra/curves/jubjub_bn183/jubjub_bn183_fields.cpp
  algebra/curves/jubjub_bn183/jubjub_bn183_g1.cpp
  algebra/curves/jubjub_bn183/jubjub_bn183_init.cpp
  algebra/curves/jubjub_bn183/jubjub_bn183_pp.cpp

  algebra/curves/jubjub_bn254/jubjub_bn254_fields.cpp
  algebra/curves/jubjub_bn254/jubjub_bn254_g1.cpp
  algebra/curves/jubjub_bn254/jubjub_bn254_init.cpp
  algebra/curves/jubjub_bn254/jubjub_bn254_pp.cpp

  algebra/curves/jubjub_ed58/jubjub_ed58_fields.cpp
  algebra/curves/jubjub_ed58/jubjub_ed58_g1.cpp
  algebra/curves/jubjub_ed58/jubjub_ed58_init.cpp
  algebra/curves/jubjub_ed58/jubjub_ed58_pp.cpp

  algebra/curves/jubjub_ed61/jubjub_ed61_fields.cpp
  algebra/curves/jubjub_ed61/jubjub_ed61_g1.cpp
  algebra/curves/jubjub_ed61/jubjub_ed61_init.cpp
  algebra/curves/jubjub_ed61/jubjub_ed61_pp.cpp

  algebra/curves/jubjub_ed97/jubjub_ed97_fields.cpp
  algebra/curves/jubjub_ed97/jubjub_ed97_g1.cpp
  algebra/curves/jubjub_ed97/jubjub_ed97_init.cpp
  algebra/curves/jubjub_ed97/jubjub_ed97_pp.cpp

  algebra/curves/jubjub_ed181/jubjub_ed181_fields.cpp
  algebra/curves/jubjub_ed181/jubjub_ed181_g1.cpp
  algebra/curves/jubjub_ed181/jubjub_ed181_init.cpp
  algebra/curves/jubjub_ed181/jubjub_ed181_pp.cpp



  algebra/curves/edwards/edwards_fields.cpp
  algebra/curves/edwards/edwards_g1.cpp
  algebra/curves/edwards/edwards_g2.cpp
  algebra/curves/edwards/edwards_init.cpp
  algebra/curves/edwards/edwards_pairing.cpp
  algebra/curves/edwards/edwards_pp.cpp

  algebra/curves/edwards58/edwards58_fields.cpp
  algebra/curves/edwards58/edwards58_g1.cpp
  algebra/curves/edwards58/edwards58_g2.cpp
  algebra/curves/edwards58/edwards58_init.cpp
  algebra/curves/edwards58/edwards58_pairing.cpp
  algebra/curves/edwards58/edwards58_pp.cpp

  algebra/curves/edwards61/edwards61_fields.cpp
  algebra/curves/edwards61/edwards61_g1.cpp
  algebra/curves/edwards61/edwards61_g2.cpp
  algebra/curves/edwards61/edwards61_init.cpp
  algebra/curves/edwards61/edwards61_pairing.cpp
  algebra/curves/edwards61/edwards61_pp.cpp

  algebra/curves/edwards97/edwards97_fields.cpp
  algebra/curves/edwards97/edwards97_g1.cpp
  algebra/curves/edwards97/edwards97_g2.cpp
  algebra/curves/edwards97/edwards97_init.cpp
  algebra/curves/edwards97/edwards97_pairing.cpp
  algebra/curves/edwards97/edwards97_pp.cpp

  algebra/curves/edwards181/edwards181_fields.cpp
  algebra/curves/edwards181/edwards181_g1.cpp
  algebra/curves/edwards181/edwards181_g2.cpp
  algebra/curves/edwards181/edwards181_init.cpp
  algebra/curves/edwards181/edwards181_pairing.cpp
  algebra/curves/edwards181/edwards181_pp.cpp


  algebra/curves/mnt/mnt4/mnt4_fields.cpp
  algebra/curves/mnt/mnt4/mnt4_g1.cpp
  algebra/curves/mnt/mnt4/mnt4_g2.cpp
  algebra/curves/mnt/mnt4/mnt4_init.cpp
  algebra/curves/mnt/mnt4/mnt4_pairing.cpp
  algebra/curves/mnt/mnt4/mnt4_pp.cpp
  algebra/curves/mnt/mnt46_common.cpp
  algebra/curves/mnt/mnt6/mnt6_fields.cpp
  algebra/curves/mnt/mnt6/mnt6_g1.cpp
  algebra/curves/mnt/mnt6/mnt6_g2.cpp
  algebra/curves/mnt/mnt6/mnt6_init.cpp
  algebra/curves/mnt/mnt6/mnt6_pairing.cpp
  algebra/curves/mnt/mnt6/mnt6_pp.cpp

  ${FF_BINARY_FIELDS}

  common/double.cpp
  common/profiling.cpp
  common/utils.cpp

  ${FF_EXTRASRCS}
)
target_link_libraries(
  ff

  ${GMPXX_LIBRARIES}
  ${GMP_LIBRARIES}

  ${PROCPS_LIBRARIES}
  ${FF_EXTRALIBS}
)
target_include_directories(
  ff
  PUBLIC ..
)

if("${BINARY_FIELDS}")
# Sodium for binary fields.
 find_path(
   SODIUM_INCLUDE_DIR
   sodium.h
   HINTS ${LIBSODIUM_INCLUDE_DIRS}
   /usr/local/include /opt/local/include /opt/include
 )
 find_library(
   SODIUM_LIBRARY
   NAMES sodium
   HINTS ${LIBSODIUM_LIBRARY_DIRS}
   /usr/local/lib /opt/local/lib /opt/lib
 )

 target_include_directories(
   ff
   PUBLIC
   ${SODIUM_INCLUDE_DIR}
 )
 target_compile_options(
   ff
   PUBLIC
   ${LIBSODIUM_CFLAGS} ${LIBSODIUM_CFLAGS_OTHER}
 )
 target_link_libraries(
   ff
   ${SODIUM_LIBRARY} ${LIBSODIUM_LDFLAGS_OTHER}
 )
endif()

install(
  DIRECTORY "" DESTINATION "include/libff"
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.tcc"
  PATTERN "tests" EXCLUDE
  PATTERN "examples" EXCLUDE
)

install(
  TARGETS ff DESTINATION lib
)

# Tests
if ("${IS_LIBFF_PARENT}")
  # Tests for fpn fields.
  add_executable(
    algebra_bilinearity_test
    EXCLUDE_FROM_ALL

    algebra/curves/tests/test_bilinearity.cpp
  )
  target_link_libraries(
    algebra_bilinearity_test

    ff
    gtest_main
  )

  add_executable(
    algebra_groups_test
    EXCLUDE_FROM_ALL

    algebra/curves/tests/test_groups.cpp
  )
  target_link_libraries(
    algebra_groups_test

    ff
    gtest_main
  )

  add_executable(
    algebra_field_utils_test
    EXCLUDE_FROM_ALL

    algebra/field_utils/tests/test_field_utils.cpp
  )
  target_link_libraries(
    algebra_field_utils_test

    ff
    gtest_main
  )

  add_executable(
    algebra_all_fields_test
    EXCLUDE_FROM_ALL

    algebra/fields/tests/test_all_fields.cpp
  )
  target_link_libraries(
    algebra_all_fields_test

    ff
    gtest_main
  )

  add_executable(
    algebra_fpn_fields_test
    EXCLUDE_FROM_ALL

    algebra/fields/tests/test_fpn_fields.cpp
  )
  target_link_libraries(
    algebra_fpn_fields_test

    ff
    gtest_main
  )

  add_executable(
          algebra_multiexp_test
          EXCLUDE_FROM_ALL

          algebra/scalar_multiplication/tests/test_multiexp.cpp
  )
  target_link_libraries(
          algebra_multiexp_test

          ff
          gtest_main
  )

  if("${BINARY_FIELDS}")
  add_executable(
    algebra_binary_fields_test
    EXCLUDE_FROM_ALL

    algebra/fields/tests/test_binary_fields.cpp
  )
  target_link_libraries(
    algebra_binary_fields_test

    ff
    gtest_main
  )
  endif()

  add_executable(
    common_test
    EXCLUDE_FROM_ALL

    common/tests/test_common.cpp
  )
  target_link_libraries(
    common_test

    ff
    gtest_main
  )

  include(CTest)
  add_test(
    NAME algebra_bilinearity_test
    COMMAND algebra_bilinearity_test
  )
  add_test(
    NAME algebra_groups_test
    COMMAND algebra_groups_test
  )
  add_test(
    NAME algebra_field_utils_test
    COMMAND algebra_field_utils_test
  )
  add_test(
    NAME algebra_all_fields_test
    COMMAND algebra_all_fields_test
  )
  add_test(
    NAME algebra_fpn_fields_test
    COMMAND algebra_fpn_fields_test
  )
  add_test(
          NAME algebra_multiexp_test
          COMMAND algebra_multiexp_test
  )
  if("${BINARY_FIELDS}")
  add_test(
    NAME algebra_binary_fields_test
    COMMAND algebra_binary_fields_test
  )
  endif()
  add_test(
    NAME common_test
    COMMAND common_test
  )

  add_dependencies(check algebra_bilinearity_test)
  add_dependencies(check algebra_groups_test)
  add_dependencies(check algebra_field_utils_test)
  add_dependencies(check algebra_all_fields_test)
  add_dependencies(check algebra_fpn_fields_test)
  add_dependencies(check algebra_multiexp_test)
  if("${BINARY_FIELDS}")
  add_dependencies(check algebra_binary_fields_test)
  endif()
  add_dependencies(check common_test)

  find_package(OpenSSL REQUIRED)
  INCLUDE_DIRECTORIES(${OPENSSL_INCLUDE_DIR})

  add_executable(
    multiexp_profile
    EXCLUDE_FROM_ALL

    algebra/scalar_multiplication/multiexp_profile.cpp
  )
  target_link_libraries(
    multiexp_profile

    ${OPENSSL_LIBRARIES}
    ff
  )

  add_executable(
          pairing_profile
          EXCLUDE_FROM_ALL

          algebra/curves/profiling/pairing_profile.cpp
  )
  target_link_libraries(
          pairing_profile

          ff
  )

  add_executable(
          addition_profile
          EXCLUDE_FROM_ALL

          algebra/curves/profiling/pairing_profile.cpp
  )
  target_link_libraries(
          addition_profile

          ff
  )

  add_dependencies(profile multiexp_profile)
  add_dependencies(profile pairing_profile)
  add_dependencies(profile addition_profile)
endif()
