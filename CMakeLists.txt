# CMake file for DISTRHO Plugins #
# ------------------------------ #

dpf_add_plugin(fildes
  TARGETS clap dssi jack lv2 vst2 vst3
  UI_TYPE cairo
  FILES_DSP
      Fildes.cpp
  FILES_UI
      Artwork.cpp
      FildesUI.cpp)

target_include_directories(
  d_cairoui PUBLIC ".")
