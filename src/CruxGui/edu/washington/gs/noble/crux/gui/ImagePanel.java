package edu.washington.gs.noble.crux.gui;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.logging.Logger;

import javax.swing.JPanel;

/**
 * The ImagePanel class is the GUI elements used to display the 
 * schematic diagram of data flow between crux components, and the
 * buttons used to display the GUI elements needed to gather the parameters
 * for each component.
 *
 * @author Charles E. Grant
 *
 */
@SuppressWarnings("serial")
class ImagePanel extends JPanel {
	private BufferedImage image;
	private static Logger logger = Logger
	.getLogger("edu.washington.gs.noble.crux.gui");

	public ImagePanel(
			CruxComponentButton createIndex,
			CruxComponentButton search,
			CruxComponentButton computeQvalues,
			CruxComponentButton percolator,
			CruxComponentButton qranker
	) {
		super();
		// read the image
		java.net.URL imgURL = getClass().getResource("/edu/washington/gs/noble/crux/gui/schematic.png");
		try {
		  image = javax.imageio.ImageIO.read(imgURL);
		}
		catch (IOException exception) {
		  logger.info(exception.getMessage());
		  System.exit(-1);
		}
		
		setSize(new Dimension(image.getWidth(), image.getHeight()));
		setMaximumSize(new Dimension(image.getWidth(), image.getHeight()));
		setMinimumSize(new Dimension(image.getWidth(), image.getHeight()));
		setPreferredSize(new Dimension(image.getWidth(), image.getHeight()));
		setLayout(null);
		createIndex.setBounds(new Rectangle(88, 61, 155, 35));
		add(createIndex);
		search.setBounds(new Rectangle(400, 61, 155, 35));
		add(search);
		computeQvalues.setBounds(new Rectangle(770, 0, 155, 35));
		add(computeQvalues);
		percolator.setBounds(new Rectangle(775, 61, 155, 35));
		add(percolator);
		qranker.setBounds(new Rectangle(770, 126, 155, 35));
		add(qranker);
	}

	public void paintComponent(Graphics g) {
		if (image == null) {
			return;
		}
		g.drawImage(image, 0, 0, getWidth(), getHeight(), this);
	}
}