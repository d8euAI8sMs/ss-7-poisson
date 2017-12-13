// poissonDlg.h : header file
//

#pragma once

#include <util/common/gui/SimulationDialog.h>

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


// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
};
