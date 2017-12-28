// poissonDlg.h : header file
//

#pragma once

#include <util/common/gui/SimulationDialog.h>
#include <util/common/plot/PlotStatic.h>

#include "model.h"

// CPoissonDlg dialog
class CPoissonDlg : public CSimulationDialog
{
// Construction
public:
	CPoissonDlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	enum { IDD = IDD_POISSON_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support

public:
    model::parameters p;
    model::relax_data d;
    model::plot_data plt;
    model::plot_data plt2;
    std::vector < std::vector < double > > P;

    virtual void OnSimulation();

// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
    double dA;
    double dE;
    double nA;
    double nE;
    BOOL A;
    BOOL E;
    bool A_bool;
    bool E_bool;
    PlotStatic plot;
    afx_msg void OnBnClickedButton1();
    afx_msg void OnBnClickedButton2();
    afx_msg void OnBnClickedCheck1();
    afx_msg void OnBnClickedCheck2();
    afx_msg void OnBnClickedCheck3();
    BOOL S;
};
