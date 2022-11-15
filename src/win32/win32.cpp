// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#if defined( _WIN32 )

#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

#include <WbemCli.h>
#include <combaseapi.h>
#pragma comment( lib, "wbemuuid.lib" )

#include <boost/assign.hpp>
#include <boost/scope_exit.hpp>
#include <boost/thread/lock_guard.hpp>
#include <boost/thread/mutex.hpp>

#include <frantic/win32/utility.hpp>

using namespace std;

namespace {
// lock used to ensure synchronized access to COM components
static boost::mutex comMutex;

// value used for initalizing COMs concurrency model
static const tagCOINIT comInitType = COINIT_MULTITHREADED;
} // namespace

namespace frantic {
namespace win32 {

namespace GetChildWindows_detail {
BOOL CALLBACK EnumChildProc( HWND hwnd, LPARAM lParam ) {
    vector<HWND>* children = reinterpret_cast<vector<HWND>*>( lParam );
    children->push_back( hwnd );
    return true;
}
} // namespace GetChildWindows_detail

frantic::tstring get_win32_computer_system_uuid::get() const {
    const boost::lock_guard<boost::mutex> comGuard( comMutex );

    /**
     * WMI or Windows Management Instrumentation allows users to manage Microsoft Windows machines. In particular it
     * allows us to query the device for various information (such as its UUID).
     *
     * COM is Microsofts component object model and is used to obtain pointers to various WMI components.
     *
     * At high level this function performs the following:
     *  Initializes COM parameters with CoInitializeEx
     *  Initializes COM process security with CoInitializeSecurity
     *  Obtains the initial locator to WMI by calling CoCreateInstance
     *  Obtains a pointer to IWbemServices (which is a called used to access WMI services
     * https://docs.microsoft.com/en-us/windows/win32/api/wbemcli/nn-wbemcli-iwbemservices and is the primary WMI
     * interface)
     *  Use the IWbemServices pointer to make a request of WMI (select uuid from Win32_ComputerSystemProduct)
     *  Get and parse the data from the WQL query.
     */

    // Initializes the COM library used for communicating with the Windows Runtime APIs
    HRESULT hRes = CoInitializeEx( nullptr, comInitType );
    if( FAILED( hRes ) ) {
        throw std::runtime_error( "verify_on_ec2_by_guid(): Unable to launch COM: " + hRes );
    }

    // Frees up any resources that are still being used and closes any remaining connections/querys
    BOOST_SCOPE_EXIT( void ) { CoUninitialize(); }
    BOOST_SCOPE_EXIT_END

    // Registers security and sets the default security values for the process
    if( ( FAILED( hRes = CoInitializeSecurity( nullptr, -1, nullptr, nullptr, RPC_C_AUTHN_LEVEL_CONNECT,
                                               RPC_C_IMP_LEVEL_IMPERSONATE, nullptr, EOAC_NONE, 0 ) ) ) ) {
        throw std::runtime_error( "verify_on_ec2_by_guid(): Unable to initialize security: " + hRes );
    }

    IWbemLocator* pLocator = nullptr;

    BOOST_SCOPE_EXIT( &pLocator ) { pLocator->Release(); }
    BOOST_SCOPE_EXIT_END

    // Creates an IWbemLocator which is an object that obtains the initial namespace pointer to the IWbemServices
    // interface for WMI on our host. Later this locator will be used to locate the pointer to the provider which
    // contains the systems UUID
    if( FAILED( hRes = CoCreateInstance( CLSID_WbemLocator, nullptr, CLSCTX_ALL, IID_PPV_ARGS( &pLocator ) ) ) ) {
        throw std::runtime_error( "verify_on_ec2_by_guid(): Unable to create a WbemLocator: " + hRes );
    }

    IWbemServices* pService = nullptr;

    BOOST_SCOPE_EXIT( &pService ) { pService->Release(); }
    BOOST_SCOPE_EXIT_END

    // connecting to the root\CIMV2 namespace
    // https://docs.microsoft.com/en-us/windows/win32/winrm/windows-remote-management-and-wmi which contains the
    // win32-computersystemproduct provider
    // https://docs.microsoft.com/en-us/windows/win32/cimwin32prov/win32-computersystemproduct
    if( FAILED( hRes = pLocator->ConnectServer( L"root\\CIMV2", nullptr, nullptr, nullptr,
                                                WBEM_FLAG_CONNECT_USE_MAX_WAIT, nullptr, nullptr, &pService ) ) ) {
        throw std::runtime_error( "verify_on_ec2_by_guid(): Unable to connect to \"CIMV2\": " + hRes );
    }

    // this is "select uuid from Win32_ComputerSystemProduct" as bytes, xored with 55, so as to not show up as literals
    // in the binary. Null terminated already.
    std::vector<wchar_t> select_uuid_xor_55 = boost::assign::list_of<wchar_t>( static_cast<wchar_t>( 68 ) )(
        static_cast<wchar_t>( 82 ) )( static_cast<wchar_t>( 91 ) )( static_cast<wchar_t>( 82 ) )( static_cast<wchar_t>(
        84 ) )( static_cast<wchar_t>( 67 ) )( static_cast<wchar_t>( 23 ) )( static_cast<wchar_t>( 66 ) )(
        static_cast<wchar_t>( 66 ) )( static_cast<wchar_t>( 94 ) )( static_cast<wchar_t>( 83 ) )( static_cast<wchar_t>(
        23 ) )( static_cast<wchar_t>( 81 ) )( static_cast<wchar_t>( 69 ) )( static_cast<wchar_t>( 88 ) )(
        static_cast<wchar_t>( 90 ) )( static_cast<wchar_t>( 23 ) )( static_cast<wchar_t>( 96 ) )( static_cast<wchar_t>(
        94 ) )( static_cast<wchar_t>( 89 ) )( static_cast<wchar_t>( 4 ) )( static_cast<wchar_t>( 5 ) )(
        static_cast<wchar_t>( 104 ) )( static_cast<wchar_t>( 116 ) )( static_cast<wchar_t>( 88 ) )(
        static_cast<wchar_t>( 90 ) )( static_cast<wchar_t>( 71 ) )( static_cast<wchar_t>( 66 ) )( static_cast<wchar_t>(
        67 ) )( static_cast<wchar_t>( 82 ) )( static_cast<wchar_t>( 69 ) )( static_cast<wchar_t>( 100 ) )(
        static_cast<wchar_t>( 78 ) )( static_cast<wchar_t>( 68 ) )( static_cast<wchar_t>( 67 ) )( static_cast<wchar_t>(
        82 ) )( static_cast<wchar_t>( 90 ) )( static_cast<wchar_t>( 103 ) )( static_cast<wchar_t>( 69 ) )(
        static_cast<wchar_t>( 88 ) )( static_cast<wchar_t>( 83 ) )( static_cast<wchar_t>( 66 ) )(
        static_cast<wchar_t>( 84 ) )( static_cast<wchar_t>( 67 ) )( static_cast<wchar_t>( 55 ) );

    for( std::vector<wchar_t>::iterator it = select_uuid_xor_55.begin(); it != select_uuid_xor_55.end(); ++it ) {
        ( *it ) ^= 55;
    }

    // This object holds a linked list of results of our "select uuid from Win32_ComputerSystemProduct" query
    IEnumWbemClassObject* pEnumerator = nullptr;

    BSTR bstr = SysAllocStringLen( select_uuid_xor_55.data(), static_cast<int>( select_uuid_xor_55.size() ) );

    BOOST_SCOPE_EXIT( &pEnumerator, &bstr ) {
        pEnumerator->Release();
        SysFreeString( bstr );
    }
    BOOST_SCOPE_EXIT_END

    // Performing the WQL (WMI Query Language Query) "select uuid from Win32_ComputerSystemProduct" against the
    // win32-computersystemproduct provider
    if( FAILED( hRes = pService->ExecQuery( L"WQL", bstr, WBEM_FLAG_FORWARD_ONLY, nullptr, &pEnumerator ) ) ) {
        throw std::runtime_error( "verify_on_ec2_by_guid(): Unable to retrive desktop monitors: " + hRes );
    }

    // Holds the current result of the query we are looking at
    IWbemClassObject* clsObj = nullptr;

    BOOST_SCOPE_EXIT( &clsObj ) { clsObj->Release(); }
    BOOST_SCOPE_EXIT_END

    // stores how many results were returned by the query. Should only be one
    int numElems;

    /**
     * enumerating over all the results from our query. The system should only return one UUID and will return the first
     * obtained UUID if more than one happen to exist.
     */
    while( ( hRes = pEnumerator->Next( WBEM_INFINITE, 1, &clsObj, (ULONG*)&numElems ) ) != WBEM_S_FALSE ) {
        if( FAILED( hRes ) )
            break;

        assert( numElems == 1 );

        VARIANT vRet;
        VariantInit( &vRet );

        BOOST_SCOPE_EXIT( &vRet ) { VariantClear( &vRet ); }
        BOOST_SCOPE_EXIT_END

        // if the result contains the property UUID return that value otherwise go to the next result
        if( SUCCEEDED( clsObj->Get( L"uuid", 0, &vRet, nullptr, nullptr ) ) && vRet.vt == VT_BSTR ) {
            return frantic::strings::to_tstring( std::wstring( vRet.bstrVal, SysStringLen( vRet.bstrVal ) ) );
        }
    }

    throw std::runtime_error( "Windows machine does not have a UUID: " + hRes );
}
} // namespace win32
} // namespace frantic

#endif // defined( _WIN32 )
