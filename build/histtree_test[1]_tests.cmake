add_test([=[HistTreeTest.DUMMY]=]  /home/mert/Hist-Tree/build/histtree_test [==[--gtest_filter=HistTreeTest.DUMMY]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[HistTreeTest.DUMMY]=]  PROPERTIES WORKING_DIRECTORY /home/mert/Hist-Tree/build SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
set(  histtree_test_TESTS HistTreeTest.DUMMY)
