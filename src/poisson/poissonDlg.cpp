// poissonDlg.cpp : implementation file
//

#include "stdafx.h"
#include "poisson.h"
#include "poissonDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CPoissonDlg dialog

using namespace plot;
using namespace util;
using namespace model;

CPoissonDlg::CPoissonDlg(CWnd* pParent /*=NULL*/)
	: CSimulationDialog(CPoissonDlg::IDD, pParent)
    , dA(1)
    , dE(1)
    , nA(10)
    , nE(10)
    , A(TRUE)
    , E(TRUE)
    , A_bool(true)
    , E_bool(true)
    , p(make_default_parameters())
    , plt(make_plot_data(RGB(100, 255, 100)))
    , plt2(make_plot_data())
    , S(FALSE)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CPoissonDlg::DoDataExchange(CDataExchange* pDX)
{
    CSimulationDialog::DoDataExchange(pDX);
    DDX_Text(pDX, IDC_EDIT12, dA);
    DDX_Text(pDX, IDC_EDIT13, dE);
    DDX_Text(pDX, IDC_EDIT14, nA);
    DDX_Text(pDX, IDC_EDIT15, nE);
    DDX_Check(pDX, IDC_CHECK1, A);
    DDX_Check(pDX, IDC_CHECK2, E);
    DDX_Text(pDX, IDC_EDIT1, p.a);
    DDX_Text(pDX, IDC_EDIT2, p.b);
    DDX_Text(pDX, IDC_EDIT3, p.k);
    DDX_Text(pDX, IDC_EDIT4, p.d);
    DDX_Text(pDX, IDC_EDIT9, p.dt);
    DDX_Text(pDX, IDC_EDIT10, p.dx);
    DDX_Text(pDX, IDC_EDIT11, p.dy);
    DDX_Text(pDX, IDC_EDIT6, p.eps);
    DDX_Text(pDX, IDC_EDIT18, p.q1);
    DDX_Text(pDX, IDC_EDIT19, p.q2);
    DDX_Text(pDX, IDC_EDIT21, p.n1);
    DDX_Text(pDX, IDC_EDIT22, p.n2);
    DDX_Text(pDX, IDC_EDIT23, p.dn1);
    DDX_Text(pDX, IDC_EDIT24, p.dn2);
    DDX_Control(pDX, IDC_PLOT, plot);
    DDX_Check(pDX, IDC_CHECK3, S);
}

BEGIN_MESSAGE_MAP(CPoissonDlg, CSimulationDialog)
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
    ON_BN_CLICKED(IDC_BUTTON1, &CPoissonDlg::OnBnClickedButton1)
    ON_BN_CLICKED(IDC_BUTTON2, &CPoissonDlg::OnBnClickedButton2)
    ON_BN_CLICKED(IDC_CHECK1, &CPoissonDlg::OnBnClickedCheck1)
    ON_BN_CLICKED(IDC_CHECK2, &CPoissonDlg::OnBnClickedCheck2)
    ON_BN_CLICKED(IDC_CHECK3, &CPoissonDlg::OnBnClickedCheck3)
END_MESSAGE_MAP()


// CPoissonDlg message handlers

BOOL CPoissonDlg::OnInitDialog()
{
	CSimulationDialog::OnInitDialog();

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	// TODO: Add extra initialization here

    adjust(p, plt);
    make_relax_data(d, p);
    P.clear();
    P.resize(d.n, std::vector < double > (d.m));

    plot.plot_layer.with(make_root_drawable(plt, {
        custom_drawable::create(make_system_painter(p, d, P)),
        plt2.plot
    }));

    plot.background = palette::brush();
    plot.triple_buffered = true;

    plot.RedrawBuffer();
    plot.SwapBuffers();
    plot.RedrawWindow();

	return TRUE;  // return TRUE  unless you set the focus to a control
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CPoissonDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CSimulationDialog::OnPaint();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CPoissonDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}


void CPoissonDlg::OnBnClickedButton1()
{
    UpdateData(TRUE);

    StartSimulationThread();
}


void CPoissonDlg::OnBnClickedButton2()
{
    StopSimulationThread();
}

void CPoissonDlg::OnSimulation()
{
    adjust(p, plt);

    make_relax_data(d, p);

    P.clear();
    P.resize(d.n, std::vector < double > (d.m));

    std::vector < plot::point < double > > hint;

    // detect field line starting points

    material_t m;
    bool border = false;
    std::vector < plot::point < size_t > > path_a, path_b;

    for (size_t l = 0; l < nE; ++l)
    {
        hint.push_back({ 1.0, (double) d.n / nE * l });
        hint.push_back({ (double) d.m - 2, (double) d.n / nE * l });
        hint.push_back({ (double) d.m / 2 - 2, (double) d.n / nE * l });
        hint.push_back({ (double) d.m / 2 + 2, (double) d.n / nE * l });
        hint.push_back({ (double) d.m / nE * l, 1.0 });
        hint.push_back({ (double) d.m / nE * l, (double) d.n - 2 });
    }

    // perform modeling

    while (m_bWorking)
    {
        relax_solve(d, p, P);
        plt.data->clear();
        plt2.data->clear();
        if (A_bool)
        {
            find_isolines(P, dA, *plt.data, d.n, d.m, p, make_material_based_stencil(d), nA);
        }
        if (E_bool)
        {
            find_field_lines(P, *plt2.data, d.n, d.m, p, make_material_based_stencil(d), hint);
        }
        plt.plot->visible = A_bool;
        plt2.plot->visible = E_bool;
        plot.RedrawBuffer();
        plot.SwapBuffers();
        Invoke([&] () { plot.RedrawWindow(); });
    }

    CSimulationDialog::OnSimulation();
}

void CPoissonDlg::OnBnClickedCheck1()
{
    UpdateData(TRUE);
    A_bool = (A != 0);
}


void CPoissonDlg::OnBnClickedCheck2()
{
    UpdateData(TRUE);
    E_bool = (E != 0);
}


void CPoissonDlg::OnBnClickedCheck3()
{
    UpdateData(TRUE);
    this->p.shift = (S != 0);
}
