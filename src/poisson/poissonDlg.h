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
    model::chasing_data d;
    model::plot_data plt;
    std::vector < std::vector < double > > P;

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
    PlotStatic plot;
};
