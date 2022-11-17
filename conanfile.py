# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import os
import typing

from conans import ConanFile, CMake

SETTINGS: dict[str, typing.Any] = {
    'os': None,
    'compiler': None,
    'build_type': None,
    'arch': ['x86_64']
}

TOOL_REQUIRES: list[str] = [
    'cmake/3.22.0',
    'thinkboxcmlibrary/1.0.0'
]

REQUIRES: list[str] = [
    'libb2/20190723',
    'bzip2/1.0.8',
    'boost/1.78.0',
    'eigen/3.4.0',
    'glog/0.5.0',
    'openexr/2.5.7',
    'tbb/2020.3',
    'tinyxml2/9.0.0',
    'utfcpp/3.2.1',
    'xxhash/0.8.1',
    'zlib/1.2.12',
    'gtest/cci.20210126',
    'xerces-c/3.2.3',
    'libe57format/2.2.0',
    'lz4/1.9.3'
]


class ThinkboxLibraryConan(ConanFile):
    name: str = 'thinkboxlibrary'
    version: str = '1.0.1'
    license: str = 'Apache-2.0'
    description: str = 'Shared code for Thinkbox Artist Tools.'
    tool_requires: list[str] = TOOL_REQUIRES
    settings: dict[str, typing.Any] = SETTINGS
    generators: str | list[str] = 'cmake_find_package'
    default_options: dict[str, typing.Any] = {
        'bzip2:build_executable': False
    }

    def requirements(self) -> None:
        for requirement in REQUIRES:
            self.requires(requirement)

        if self.settings.os == 'Linux':
            self.requires('icu/71.1')

    def build(self) -> None:
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
    
    def imports(self) -> None:
        # Copy tbb DLLs to the UnitTests binary directory
        self.copy('*.dll', dst='UnitTests/Release', src='bin')

    def export_sources(self) -> None:
        self.copy('**.h', src='', dst='')
        self.copy('**.hpp', src='', dst='')
        self.copy('**.cpp', src='', dst='')
        self.copy('**.cc', src='', dst='')
        self.copy('UnitTests/TestInputs/*', src='', dst='')
        self.copy('UnitTests/CMakeLists.txt', src='', dst='')
        self.copy('CMakeLists.txt', src='', dst='')
        self.copy('NOTICE.txt', src='', dst='')
        self.copy('LICENSE.txt', src='', dst='')
        self.copy('THIRD-PARTY-LICENSES', src='', dst='')

    def package(self) -> None:
        cmake = CMake(self)
        cmake.install()

        with open(os.path.join(self.source_folder, 'NOTICE.txt'), 'r', encoding='utf8') as notice_file:
            notice_contents = notice_file.readlines()
        with open(os.path.join(self.source_folder, 'LICENSE.txt'), 'r', encoding='utf8') as license_file:
            license_contents = license_file.readlines()
        with open(os.path.join(self.source_folder, 'THIRD-PARTY-LICENSES'), 'r', encoding='utf8') as third_party_file:
            third_party_contents = third_party_file.readlines()
        os.makedirs(os.path.join(self.package_folder, 'licenses'), exist_ok=True)
        with open(os.path.join(self.package_folder, 'licenses', 'LICENSE'), 'w', encoding='utf8') as cat_license_file:
            cat_license_file.writelines(notice_contents)
            cat_license_file.writelines(license_contents)
            cat_license_file.writelines(third_party_contents)

    def deploy(self) -> None:
        self.copy("*", dst="lib", src="lib")
        self.copy("*", dst="include", src="include")

    def package_info(self) -> None:
        self.cpp_info.libs = ["thinkboxlibrary"]
